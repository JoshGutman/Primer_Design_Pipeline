import subprocess
import argparse
import glob
import os
from Bio import SeqIO

from get_seqs import get_seqs
from get_genomes import get_genomes
from generate_primers import generate_primers


#TODO
# Output all mis-hits and non-target hits to file
# One file for each primer mis-hit and each primer non-target-hit
# E.g. 244_reverse_mis_hits.txt, 244_reverse_non_target_hits.txt


def primer_design_pipeline(target, directory, config_file, target_list, lower, upper, ignore):

    combine_seqs(directory)
    
    seqs, mis_hits, non_target_hits = get_seqs(config_file, target, directory, target_list, lower, upper)
    genomes = get_genomes(target, directory)
    primers = generate_primers(seqs, genomes, ignore)

    combos = get_combos(primers, lower, upper)
    output_candidate_primers(combos, primers, mis_hits, non_target_hits)

    print(primers)





# Create fasta used to make blast database
def combine_seqs(directory):

    '''
    if os.path.isfile("combined.seqs"):
        return
    '''
    
    with open("combined.seqs", "w") as f:

        for file in glob.glob(os.path.join(directory, "*.fasta")):
            f.write(">" + os.path.basename(file).replace(".fasta", "") + "\n")

            # BioPython dependancy
            with open(file) as temp_file:
                for record in SeqIO.parse(temp_file, "fasta"):
                    f.write(str(record.seq) + "\n")

    subprocess.run("makeblastdb -in combined.seqs -dbtype nucl > /dev/null 2>&1", shell=True)



def get_combos(primers, lower, upper):

    out = {}
    all_combos = {}
    all_combos["no matches"] = None
    forwards = []
    reverses = []

    for primer in primers:
        if "forward" in primer:
            forwards.append(primer)
            all_combos[primer] = []
        else:
            reverses.append(primer)

    for f in forwards:
        f_val = int(f.split("_")[0])
        
        for r in reverses:
            r_val = int(r.split("_")[0])

            combo_range = abs(f_val - r_val)

            all_combos[f].append((r, combo_range))
            
            if lower <= combo_range <= upper:

                if f in out:
                    out[f].append((r, combo_range))
                else:
                    out[f] = [(r, combo_range)]


    if len(out) != 0:
        print(out)
        print()
        return out

    else:
        return all_combos





def output_candidate_primers(combos, primers, mis_hits, non_target_hits):

    with open("candidate_primers.txt", "w") as outfile:
        #outfile.write("Forward name\tReverse name\tMis-hits sum\tNon-target hits sum\t# forward degens\t # reverse degens\tForward sequence\tReverse Sequence\n")
        outfile.write("Forward name\tReverse name\tMax mis-hit (ID, Length)\tMax non-target-hit (ID, length)\t# forward degens\t # reverse degens\tForward sequence\tReverse Sequence\n")

        if "no_matches" in combos:
            # No combos were found in the specified range
            #TODO
            pass

        else:
            for forward in combos:

                for data in combos[forward]:
                    reverse = data[0]
                    
                    vals = []
                    
                    vals.append(forward)
                    vals.append(reverse)
                    #vals.append(len(mis_hits[forward]) + len(mis_hits[reverse]))
                    #vals.append(len(non_target_hits[forward]) + len(non_target_hits[reverse]))
                    vals.append(mis_hits[0])
                    vals.append(non_target_hits[0])
                    vals.append(get_number_degens(primers[forward]))
                    vals.append(get_number_degens(primers[reverse]))
                    vals.append(primers[forward])
                    vals.append(primers[reverse])

                    outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(*vals))
        




def get_number_degens(sequence):

    out = 0
    for base in sequence.upper():
        if base not in "ACGT":
           out += 1 

    return out


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="primer_design_pipeline.py")
    parser.add_argument("-t", "--target", help="[REQUIRED] Path to target multifasta", required=True)
    parser.add_argument("-d", "--directory", help="[REQUIRED] Path to reference fastas", required=True)
    parser.add_argument("-c", "--config", help="[REQUIRED] Path to primer3 config file", required=True)
    parser.add_argument("-g", "--genomes", help="[REQUIRED] Path to .txt file with target genomes", required=True)
    parser.add_argument("-l", "--lower", help="Lower bound of amplicon size", type=int, default=150)
    parser.add_argument("-u", "--upper", help="Upper bound of amplicon size", type=int, default=250)
    parser.add_argument("-i", "--ignore", help="Threshold percentage to consider degens", type=int, default=95)

    args = parser.parse_args()

    primer_design_pipeline(args.target, args.directory, args.config, args.genomes, args.lower, args.upper, args.ignore)
