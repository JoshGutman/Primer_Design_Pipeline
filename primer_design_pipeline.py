import argparse

from get_seqs import get_seqs
from get_genomes import get_genomes
from generate_primers import generate_primers


def primer_design_pipeline(target, directory, config_file, target_list, lower, upper):

    combine_seqs(directory)
    seqs = get_seqs(config_file, target, directory, target_list, lower, upper)
    genomes = get_genomes(target, directory)
    primers = generate_primers(seqs, genomes)

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



if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="primer_design_pipeline.py")
    parser.add_argument("-t", "--target", help="[REQUIRED] Path to target multifasta", required=True)
    parser.add_argument("-d", "--directory", help="[REQUIRED] Path to reference fastas", required=True)
    parser.add_argument("-c", "--config", help="[REQUIRED] Path to primer3 config file", required=True)
    parser.add_argument("-g", "--genomes", help="[REQUIRED] Path to .txt file with target genomes", required=True)
    parser.add_argument("-l", "--lower", help="Lower bound of amplicon size", type=int, default=150)
    parser.add_argument("-u", "--upper", help="Upper bound of amplicon size", type=int, default=250)

    args = parser.parse_args()

    primer_design_pipeline(args.target, args.directory, args.config, args.genomes, args.lower, args.upper)
