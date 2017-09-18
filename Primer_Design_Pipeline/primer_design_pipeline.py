import subprocess
import argparse
import pickle
import sys
import os

from .setup import Constants, init
from .get_tm import get_tm
from .get_primers import get_primers
from .get_genomes import get_genomes
from .get_degens import get_degens
from .combo_funcs import *
from .find_primer_conflicts import find_primer_conflicts, blast_all_primers


# Driver
def primer_design_pipeline(target_file, directory, config_file, target_list,
                           reference_fasta, lower, upper, ignore, oligo_conc,
                           na_conc, mg_conc, project_dir, keep):

    target_list = os.path.abspath(target_list)
    reference_fasta = os.path.abspath(reference_fasta)
    config_file = os.path.abspath(config_file)

    targets = init(target_file, directory, config_file, target_list,
                   reference_fasta, project_dir)

    all_combos = []
    best_combos = []
    new_dirs = []
    for target in targets:

        dir_name = os.path.splitext(target)[0]
        os.mkdir(dir_name)
        new_dirs.append(dir_name)
        
        subprocess.run("cp {} {}".format(target, dir_name), shell=True)
        subprocess.run("cp {} {}".format(Constants.config_file, dir_name),
                       shell=True)
        os.chdir(dir_name)
        config_file = os.path.abspath(os.path.basename(config_file))


        primers, mis_hits, non_target_hits = get_primers(config_file, target, directory, lower, upper)

        genomes = get_genomes(target)
        get_degens(primers, genomes, ignore)

        blast_all_primers("alignment_blast_in.fasta")
        find_primer_conflicts("alignment_blast_in.fasta")

        combos = get_combos(primers, lower, upper)

        for combo in combos:
            combo.forward.tm = get_tm(combo.forward.sequence,
                                      oligo_conc, na_conc, mg_conc)
            combo.reverse.tm = get_tm(combo.reverse.sequence,
                                      oligo_conc, na_conc, mg_conc)

            combo.forward.set_ordering_info(target, combo.amplicon)
            combo.reverse.set_ordering_info(target, combo.amplicon)

            combo.target = target
            all_combos.append(combo)

        score_combos(primers, combos)

        output_combos(combos, "candidate_primers.txt")
        #output_ordering_info(combos)
        best_combos.append(choose_best_combos(combos))

        os.chdir("..")

    for combo in best_combos:
        output_combos(combo, "best_primers.txt")

    pickle_combos(all_combos, fix_imports=False)

    if not keep:
        remove_excess_files(new_dirs)


def remove_excess_files(directories):

    # Get rid of files in sub-directories
    for directory in directories:
        os.chdir(directory)
        
        subprocess.run("rm {}*".format(FileNames.conflict_blast_input), shell=True)
        subprocess.run("rm {}*".format(directory), shell=True)

        os.remove(FileNames.neben_output)
        os.remove(FileNames.modified_config_file)
        os.remove(FileNames.primer3_output)
        os.remove(FileNames.target_blast)
        os.remove(FileNames.non_target_blast)
        os.remove(os.path.basename(Constants.config_file))
        os.remove(FileNames.muscle_input)
        os.remove(FileNames.muscle_output)

        os.chdir("..")

    # Get rid of all .fasta files from multifasta in outer-most directory
    for directory in directories:
        os.remove(directory + ".fasta")

    subprocess.run("rm {}*".format(Constants.combined_seqs), shell=True)
    subprocess.run("rm {}*".format(Constants.target_db), shell=True)
    subprocess.run("rm {}*".format(Constants.non_target_db), shell=True)
        

def pickle_combos(all_combos):
    with open(FileNames.pickled_combos, "wb") as outfile:
        pickle.dump(all_combos, outfile)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="primer_design_pipeline.py")
    parser.add_argument("-t", "--target", help="[REQUIRED] Path to target multifasta", required=True)
    parser.add_argument("-d", "--directory", help="[REQUIRED] Path to reference fastas", required=True)
    parser.add_argument("-c", "--config", help="[REQUIRED] Path to primer3 config file", required=True)
    parser.add_argument("-g", "--genomes", help="[REQUIRED] Path to .txt file with target genomes", required=True)
    parser.add_argument("-r", "--reference", help="Path to amplicon reference assembly .fasta", default="combined.seqs")
    parser.add_argument("-l", "--lower", help="Lower bound of amplicon size", type=int, default=150)
    parser.add_argument("-u", "--upper", help="Upper bound of amplicon size", type=int, default=250)
    parser.add_argument("-i", "--ignore", help="Threshold percentage to consider degens", type=float, default=98)
    parser.add_argument("-ol", "--oligo_conc", help="Oligo concentration (μM) for calculating Tm", type=float, default=.25)
    parser.add_argument("-na", "--na_conc", help="Na+ concentration (mM) for calculating Tm", type=float, default=50)
    parser.add_argument("-mg", "--mg_conc", help="Mg++ concentration (mM) for calculating Tm", type=float, default=0)
    parser.add_argument("-k", "--keep", help="Keep all temporary files", type=bool, default=False) #TODO

    project_dir = os.path.dirname(os.path.abspath(sys.argv[0]))

    args = parser.parse_args()
    primer_design_pipeline(args.target, args.directory, args.config, args.genomes, args.reference, args.lower, args.upper, args.ignore, args.oligo_conc, args.na_conc, args.mg_conc, project_dir, args.keep)