import argparse

'''
from get_seqs import get_seqs
from get_genomes import get_genomes
from generate_primers import generate_primers


def primer_design_pipeline(target, directory, config_file, lower, upper):

    seqs = get_seqs(config_file, target, lower, upper)
    genomes = get_genomes(target, directory)
    primers = generate_primers(seqs, genomes)

    print(primers)
'''


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="primer_design_pipeline.py")
    parser.add_argument("-t", "--target", help="[REQUIRED] Path to target multifasta", required=True)
    parser.add_argument("-d", "--directory", help="[REQUIRED] Path to reference fastas", required=True)
    parser.add_argument("-c", "--config", help="[REQUIRED] Path to primer3 config file", required=True)
    parser.add_argument("-l", "--lower", help="Lower bound of amplicon size", type=int, default=150)
    parser.add_argument("-u", "--upper", help="Upper bound of amplicon size", type=int, default=250)

    args = parser.parse_args()

    primer_design_pipeline(args.target, args.directory, args.config, args.lower, args.upper)
