"""Uses BLAST to find conflicts between candidate primers.

Also includes function blast_all_primers to BLAST candidate primers
against a database of all reference genomes.

"""
import os
import itertools
import subprocess

from Bio import Seq, SeqIO

from .setup import Constants, FileNames

def find_primer_conflicts(primer_fasta):
    """Finds conflicts between all candidate primers

    Makes a blast database from all of the primer candidates, and then BLASTs
    the candidates against that database. Keeps track of all conflics of
    lengths in the range [0, 10].

    Args:
        primer_fasta (str): Path to .fasta file with candidate primers.

    Returns:
        dict: Dictionary of conflicts:

            key (int): alignment length of conflict.
            value (:obj:`list` of :obj:`list` of :obj:`str`): List of each
            conflict, where each conflict is a list where item 0 conflicts
            with item 1.

    Notes:
        Doesn't account for the fact that values of higher alignment lengths
        should be encapsulated in lower alignment lengths. For example,
        conflicts[8] should contain the conflicts that occur at alignment
        lengths of 8, 9, and 10, but only contain the the conflicts
        that occur at length 8.

    """
    primer_path = os.path.abspath(primer_fasta)

    all_seqs = _get_all_sequences(primer_path)

    with open(primer_path, "w") as outfile:
        for key in all_seqs:
            for seq in all_seqs[key]:
                outfile.write(">{}\n{}\n".format(key, seq))
    
    subprocess.run("makeblastdb -in {} -dbtype nucl > /dev/null 2>&1".format(primer_path), shell=True)
    subprocess.run("blastn -task blastn -query {} -db {} -perc_identity 100 -dust no -evalue 20 -word_size 4 -outfmt 6 -out primer_conflicts_blast.out > /dev/null 2>&1".format(primer_path, primer_path), shell=True)

    out = _parse_blast_output("primer_conflicts_blast.out")
    #os.remove("primer_conflicts_blast.out")

    #output_conflicts(out)

    return out


def _parse_blast_output(output):
    """Gets information from BLAST output

    Gets conflicting primers pairs at alignment lengths of [0, 10].

    Args:
        output (str): Path to BLAST output file.

    Returns:
        dict: Dictionary of conflicts:

            key (int): alignment length of conflict.
            value (:obj:`list` of :obj:`list` of :obj:`str`): List of each
            conflict, where each conflict is a list where item 0 conflicts
            with item 1.

    Notes:
        Doesn't account for the fact that values of higher alignment lengths
        should be encapsulated in lower alignment lengths. For example,
        conflicts[8] should contain the conflicts that occur at alignment
        lengths of 8, 9, and 10, but in reality, only contain the the conflicts
        that occur at length 8.

    """
    conflicts = {}

    # If there are no conflicts at a certain length, conflicts[length] == []
    for i in range(1, 11):
        conflicts[i] = []

    duplicates = set()

    with open(output, "rU") as f:

        for line in f:
            fields = line.strip().split()

            def _conflict_is_unique():
                # fields[3]=length, fields[0]=primer #1, fields[1]=primer #2
                # What is fields[5]?
                if int(fields[3]) > 10 or int(fields[5]) != 0:
                    return False
                if (fields[0], fields[1]) in duplicates:
                    return False
                if (fields[1], fields[0]) in duplicates:
                    return False
                if fields[0] == fields[1]:
                    return False
                return True

            if _conflict_is_unique():
                conflicts[int(fields[3])].append([fields[0], fields[1]])
                duplicates.add((fields[0], fields[1]))

    return conflicts


def output_conflicts(conflicts):
    """Outputs conflicts to a file.

    Writes conflicting primer pairs to "primer_conflicts.txt".

    Args:
        conflicts (dict): Dictionary of conflicts.

    Returns:
        None

        Creates file "primer_conflicts.txt".

    Notes:
        `conflicts` is outputted from _parse_blast_output.

    """
    with open("primer_conflicts.txt", "w") as out:
        for length in conflicts:
            if len(conflicts[length]) > 0:
                out.write("Problematic primers with"
                          "an alignment length of {}:\n".format(length))
                for i in range(10, length-1, -1):
                    if len(conflicts[i]) > 0:
                        for data in conflicts[i]:
                            out.write("{}\t{}\n".format(data[0], data[1]))
                out.write("\n\n")


def _get_all_sequences(fasta):
    to_write = {}
    with open(fasta) as infile:
        for record in SeqIO.parse(infile, "fasta"):
            to_write[str(record.id)] = _expand_degenerate_sequence(
                str(record.seq))
    return to_write


def _expand_degenerate_sequence(seq):
    d = Seq.IUPAC.IUPACData.ambiguous_dna_values
    return list(map("".join, itertools.product(*map(d.get, seq))))
                

def blast_all_primers(primer_fasta):
    """BLAST all candidate primers against database of all reference genomes.

    Args:
        primer_fasta (str): Path to .fasta file with candidate primers.

    Returns:
        None

        Creates file "primers.blast.out".

    """
    subprocess.run("blastall -p blastn -i {} -d {} -m 8 -e 10 -o primers.blast.out".format(primer_fasta, Constants.combined_seqs), shell=True)
