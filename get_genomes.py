import glob
import os
import subprocess
from Bio import SeqIO

"""
Will tblastn ever be run?
Is Muscle needed? (Probably not)
Is genome length of 500 a magic number? (Line 89)
"""


# Driver
# Takes in target.fasta and directory with reference genomes
def get_genomes(target, directory):

    directory = os.path.abspath(directory)

    combine_seqs(directory)
        
    run_blast(target, "blastn")
    parse_blast_output(target)

    subprocess.run("muscle -in tmp_muscle_in -out tmp_muscle_out.fasta", shell=True)

    out = parse_muscle_output()


    # Remove excess files
        # Remove muscle files
    for file in glob.glob("tmp_muscle_*"):
        os.remove(file)

        # Remove blast database files
    for file in glob.glob("combined.seqs*"):
        os.remove(file)

        # Remove target fasta variants (leave original target.fasta)
    os.remove("{}.blast.out".format(target))
    os.remove("{}.for".format(target))
    os.remove("{}.rev".format(target))
    
    
    

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



# Run blast
def run_blast(target, blast_type):

    if blast_type == "blastn":
        reduced_name = os.path.basename(target).replace(".fasta", "")
        subprocess.run('blastn -task blastn -query {} -db combined.seqs -out {}.blast.out -dust no -num_alignments 20000 -outfmt "7 std sseq"'.format(target, reduced_name), shell=True)
    # Is tblastn needed?
    else:
        reduced_name = os.path.basename(target).replace(".pep", "")
        subprocess.run('tblastn -query {} -db combined.seqs -out {}.blast.out -seg no -num_alignments 20000 -outfmt "7 std sseq"'.format(target, reduced_name), shell=True)



# Get seqs outputted from blast
def parse_blast_output(target):

    with open("{}.blast.out".format(os.path.basename(target).replace(".fasta", "")), "rU") as infile, open("tmp_muscle_in", "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                continue
            else:
                fields = line.split()
                # Is 500 a magic number?
                if len(fields[12]) == 500:
                    outfile.write(">" + fields[1] + "\n")
                    outfile.write(fields[12] + "\n")


# Is muscle needed?
def parse_muscle_output():

    out = []

    with open("tmp_muscle_out.fasta") as f:
        for record in SeqIO.parse(f, "fasta"):
            out.append(str(record.seq))
                    

    return out



if __name__ == "__main__":

    import sys

    args = sys.argv

    print(get_genomes(args[1], args[2]))
