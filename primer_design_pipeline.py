import subprocess
import argparse
import glob
import sys
import os
from Bio import SeqIO

from get_tm import get_tm
from get_primers import get_primers
from get_genomes import get_genomes
from get_degens import get_degens
from find_primer_conflicts import find_primer_conflicts, blast_all_primers


# Driver
def primer_design_pipeline(target, directory, config_file, target_list, reference_fasta, lower, upper, ignore, oligo_conc, na_conc, mg_conc, project_dir):

    combine_seqs(directory)
    
    primers, mis_hits, non_target_hits = get_primers(config_file, target, directory, target_list, lower, upper)
    genomes = get_genomes(target, directory)
    get_degens(primers, genomes, ignore)

    blast_all_primers("alignment_blast_in.fasta", "combined.seqs")
    find_primer_conflicts("alignment_blast_in.fasta")

    combos = get_combos(primers, lower, upper, reference_fasta, project_dir)
    
    output_candidate_primers(combos, target, [oligo_conc, na_conc, mg_conc])

    print(primers)


# Create fasta used to make blast database
def combine_seqs(directory):

    if os.path.isfile("combined.seqs"):
        return

    with open("combined.seqs", "w") as f:

        for file in glob.glob(os.path.join(directory, "*.fasta")):
            f.write(">" + os.path.basename(file).replace(".fasta", "") + "\n")

            # BioPython dependancy
            with open(file) as temp_file:
                for record in SeqIO.parse(temp_file, "fasta"):
                    f.write(str(record.seq) + "\n")

    subprocess.run("makeblastdb -in combined.seqs -dbtype nucl > /dev/null 2>&1", shell=True)



def get_combos(primers, lower, upper, reference_fasta, project_dir):
    out = set()

    forwards = []
    reverses = []

    for primer in primers:
        if primer.orientation == "forward":
            forwards.append(primer)
        elif primer.orientation == "reverse":
            reverses.append(primer)

    for f in forwards:
        for r in reverses:
            combo_range = abs(f.value - r.value)
            if lower <= combo_range <= upper:
                temp_combo = Combo(f, r)
                temp_combo.set_amplicon(reference_fasta, project_dir)
                out.add(temp_combo)

    return out




def output_candidate_primers(combos, target, tm_args):

    with open("candidate_primers.txt", "w") as outfile:

        if len(combos) == 0:
            # No combos were found in the specified range
            #TODO
            pass

        else:
            for combo in combos:
                combo.forward.tm = get_tm(combos.forward.sequence, *tm_args)
                combo.reverse.tm = get_tm(combos.reverse.sequence, *tm_args)

                outfile.write(combo.name + "\n")
                outfile.write("--------------------------------------------------------------------------------------\n")
                
                # Name, mis-hit, non-target hit, degens, sequence, tm
                outfile.write("Name\t\tMax mis-hit\t\tMax non-target hit\t\t# degens\t\tsequence\t\t[Min,Max,Avg] Tm\n")
                outfile.write("{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\n".format(combo.forward.name,
                                                                          combo.forward.max_mis_hit,
                                                                          combo.forward.max_non_target_hit,
                                                                          combo.forward.num_degens,
                                                                          combo.forward.sequence,
                                                                          combo.forward.tm))

                outfile.write("{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\n".format(combo.reverse.name,
                                                                          combo.reverse.max_mis_hit,
                                                                          combo.reverse.max_non_target_hit,
                                                                          combo.reverse.num_degens,
                                                                          combo.reverse.sequence,
                                                                          combo.reverse.tm))

                combo.forward.set_ordering_info(target, combo.amplicon)
                combo.reverse.set_ordering_info(target, combo.amplicon)

                outfile.write("\nOrdering information:\n")
                # Target, Primer, Combined_Name, Primer (5'-3'), final_name, UT + Sequnece, To order, Tm, Amplicon, Amplicon Length, Amplicon length + UT
                outfile.write(combo.forward.ordering_info + "\n")
                outfile.write(combo.reverse.ordering_info + "\n")

                # Target, Amplicon
                outfile.write("{},{}\n".format(target, combo.amplicon))

                outfile.write("\n\n\n")




class Combo:

    def __init__(self, forward_primer, reverse_primer):
        self.forward = forward_primer
        self.reverse = reverse_primer
        self.name = "{} - {}".format(self.forward.name, self.reverse.name)
        self.amplicon = None
        self.target = None
        self.primer_name = None
        self.combined_name = None
        self.sequence = None

    def set_amplicon(self, reference_fasta, project_dir):
        
        subprocess.run("{}/neben_linux_64 --primers {}:{} {} > amplicon".format(project_dir,
                                                                                self.forward.sequence,
                                                                                self.reverse.sequence,
                                                                                reference_fasta),shell=True)
        with open("amplicon", "rU") as f:
            out = f.readline()
            
        if out == "":
            self.amplicon = "None found"
        else:
            self.amplicon = out.split()[3]



if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="primer_design_pipeline.py")
    parser.add_argument("-t", "--target", help="[REQUIRED] Path to target multifasta", required=True)
    parser.add_argument("-d", "--directory", help="[REQUIRED] Path to reference fastas", required=True)
    parser.add_argument("-c", "--config", help="[REQUIRED] Path to primer3 config file", required=True)
    parser.add_argument("-g", "--genomes", help="[REQUIRED] Path to .txt file with target genomes", required=True)
    parser.add_argument("-r", "--reference", help="[REQUIRED] Path to amplicon reference genome .fasta", required=True)
    parser.add_argument("-l", "--lower", help="Lower bound of amplicon size", type=int, default=150)
    parser.add_argument("-u", "--upper", help="Upper bound of amplicon size", type=int, default=250)
    parser.add_argument("-i", "--ignore", help="Threshold percentage to consider degens", type=int, default=95)
    parser.add_argument("-ol", "--oligo_conc", help="Oligo concentration (μM) for calculating Tm", type=float, default=.25)
    parser.add_argument("-na", "--na_conc", help="Na+ concentration (mM) for calculating Tm", type=float, default=50)
    parser.add_argument("-mg", "--mg_conc", help="Mg++ concentration (mM) for calculating Tm", type=float, default=0)
    parser.add_argument("-k", "--keep", help="Keep all temporary files", type=bool, default=False) #TODO

    project_dir = os.path.dirname(os.path.abspath(sys.argv[0]))

    args = parser.parse_args()
    primer_design_pipeline(args.target, args.directory, args.config, args.genomes, args.reference, args.lower, args.upper, args.ignore, args.oligo_conc, args.na_conc, args.mg_conc, project_dir)
