import glob
import os
import subprocess
from Bio import SeqIO

"""
Will tblastn ever be run?
Is genome length of 500 a magic number? (Line 89)
"""


# Driver
# Takes in target.fasta and directory with reference genomes
def get_genomes(target, directory):

    directory = os.path.abspath(directory)
        
    run_blast(target, "blastn")
    parse_blast_output(target)

    subprocess.run("muscle -in tmp_muscle_in -out tmp_muscle_out.fasta > /dev/null 2>&1", shell=True)

    out = parse_muscle_output()

    '''
    # Remove excess files
        # Remove muscle files
    for file in glob.glob("tmp_muscle_*"):
        os.remove(file)

        # Remove blast database files
    for file in glob.glob("combined.seqs*"):
        os.remove(file)
    

        # Remove target fasta variants (leave original target.fasta)
    target_name = os.path.splitext(target)[0]
    os.remove("{}.blast.out".format(target_name))
    os.remove("{}.for".format(target_name))
    os.remove("{}.rev".format(target_name))
    '''

    return out
    



# Run blast
def run_blast(target, blast_type):

    if blast_type == "blastn":
        reduced_name = os.path.basename(target).replace(".fasta", "")
        subprocess.run('blastn -task blastn -query {} -db {} -out {}.blast.out -dust no -num_alignments 20000 -outfmt "7 std sseq"'.format(target, COMBINED_SEQS, reduced_name), shell=True)
    # Is tblastn needed?
    else:
        reduced_name = os.path.basename(target).replace(".pep", "")
        subprocess.run('tblastn -query {} -db {} -out {}.blast.out -seg no -num_alignments 20000 -outfmt "7 std sseq"'.format(target, COMBINED_SEQS, reduced_name), shell=True)


def parse_blast_output(target):

    data = []

    with open("{}.blast.out".format(os.path.basename(target).replace(".fasta", "")), "rU") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            else:
                fields = line.split()
                data.append((fields[1], fields[12]))

    data.sort(key=lambda tup: len(tup[1]), reverse=True)
    largest = len(data[0][1])

    with open("tmp_muscle_in", "w") as outfile:
        for item in data:
            if (len(item[1]) / largest) > .9:
                outfile.write(">{}\n".format(item[0]))
                outfile.write("{}\n".format(item[1]))
        

def parse_muscle_output():

    out = []

    with open("tmp_muscle_out.fasta") as f:
        for record in SeqIO.parse(f, "fasta"):
            out.append(str(record.seq))
                    
    return out
