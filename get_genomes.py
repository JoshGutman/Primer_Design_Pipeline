"""Gets the sequences that will be used to look for degeneracies.

The driver function, get_genomes, uses BLAST, sort, and MUSCLE to get the
sequences.

"""
import os
import subprocess
from Bio import SeqIO

from setup import Constants, FileNames


def get_genomes(target):
    """Get the genomes that will be used to look for degeneracies.

    BLAST the target sequence against the combined.seqs BLAST database. Then,
    get the unique subject IDs and sequences from the BLAST output. Finally,
    use MUSCLE to align all of the sequences.

    Args:
        target (str): Path to target .fasta file

    Returns:
        :obj:`list` of :obj:`str`: Aligned sequences

    Notes:
        Driver function

    """
    # "/path/to/target.fasta" -> "target"
    reduced_name = os.path.splitext(os.path.basename(target))[0]

    # BLAST target against BLAST database and get output
    run_blast(target, reduced_name)
    blast_data = parse_blast_output(reduced_name)

    # Create MUSCLE input, run MUSCLE, and get results
    create_muscle_input(blast_data)
    subprocess.run("muscle -in {} -out {} > /dev/null 2>&1"
                   .format(FileNames.muscle_input, FileNames.muscle_output),
                   shell=True)
    return parse_muscle_output()


def run_blast(target, reduced_name):
    """BLAST target file against BLAST database and get rid of duplciates.

    Args:
        target (str): Path to target .fasta file.
        reduced_name (str): target basename without extension.

    Returns:
        None.

        Outputs file names <reduced_name>.blast.out

    Notes:
        Linux sort command is used to both sort the output file in place and get
        rid of any duplicates (using the -u flag).

    """
    subprocess.run('blastn -task blastn -query {} -db {} -out {}.blast.out '
                   '-dust no -num_alignments 20000 -outfmt "7 std sseq"'
                   .format(target, Constants.combined_seqs, reduced_name),
                   shell=True)

    # Sort output in place and get rid of duplicates
    subprocess.run("sort -u -k 2,2 {}.blast.out -o {}.blast.out"
                   .format(reduced_name, reduced_name), shell=True)


def parse_blast_output(reduced_name):
    """Get every subject ID and sequence from BLAST output.

    Args:
        reduced_name (str): Name of the target .fasta file without extensions.

    Returns:
        :obj:`list` of :obj:`tuple` of :obj:`str`: List of tuples, each tuple
        being a subject ID and subject sequence from blast output.

    """
    data = []
    with open("{}.blast.out".format(reduced_name), "rU") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            else:
                fields = line.split()
                data.append((fields[1], fields[12]))
    return data


def create_muscle_input(blast_data):
    """Creates an input file for MUSCLE based on BLAST output data.

    Args:
        blast_data (:obj:`list` of :obj:`tuple` of :obj:`str`): Output from
            parse_blast_output.

    Returns:
        None.

        Creates input file for MUSCLE.

    Notes:
        Input file name defined in setup.py in class FileNames.

    """
    # Sort blast data based on the length of sequences
    blast_data.sort(key=lambda tup: len(tup[1]), reverse=True)
    largest_len = len(blast_data[0][1])

    with open(FileNames.muscle_input, "w") as outfile:
        for record in blast_data:
            # If the sequence is at least 90% the size of the largest sequence
            if (len(record[1]) / largest_len) > .9:
                outfile.write(">{}\n".format(record[0]))
                outfile.write("{}\n".format(record[1]))


def parse_muscle_output():
    """Get genomes from MUSCLE output.

    Parses the MUSCLE output file using Bio.SeqIO to get sequences.

    Args:
        None.

    Returns:
        :obj:`list` of :obj:`str`: List of genomes.

    Notes:
        MUSCLE output file name degined in setup.py in class FileNames.

    """
    out = []
    with open(FileNames.muscle_output) as infile:
        for record in SeqIO.parse(infile, "fasta"):
            out.append(str(record.seq))
    return out
