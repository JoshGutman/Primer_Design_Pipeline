import glob
import os
import subprocess
from Bio import SeqIO

"""
Will tblastn ever be run?
"""

# Driver
def extract_seqs(target, path):

    path = os.path.abspath(path)

    if not os.path.isfile("combined.seqs"):
        combine_seqs(path)
        subprocess.run("makeblastdb -in combined.seqs -dbtype nucl > /dev/null 2>&1", shell=True)

    run_blast(target, "blastn")

    return parse_blast_output(target)

    

# Create fasta used to make blast database
def combine_seqs(path):
    
    with open("combined.seqs", "w") as f:

        for file in glob.glob(os.path.join(path, "*.fasta")):
            f.write(">" + os.path.basename(file).replace(".fasta", "") + "\n")
            
            with open(file) as temp_file:
                for record in SeqIO.parse(temp_file, "fasta"):
                    f.write(str(record.seq) + "\n")



# Run blast
def run_blast(target, blast_type):

    if blast_type == "blastn":
        reduced_name = os.path.basename(target).replace(".fasta", "")
        subprocess.run('blastn -task blastn -query {} -db combined.seqs -out {}.blast.out -dust no -num_alignments 20000 -outfmt "7 std sseq"'.format(target, reduced_name), shell=True)
    else:
        reduced_name = os.path.basename(target).replace(".pep", "")
        subprocess.run('tblastn -query {} -db combined.seqs -out {}.blast.out -seg no -num_alignments 20000 -outfmt "7 std sseq"'.format(target, reduced_name), shell=True)



# Get seqs outputted from blast
def parse_blast_output(target):

    out = []

    with open("{}.blast.out".format(os.path.basename(target).replace(".fasta", "")), "rU") as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                out.append(line.split()[12])

    return out


if __name__ == "__main__":

    import sys

    args = sys.argv

    print(extract_seqs(args[1], args[2]))
