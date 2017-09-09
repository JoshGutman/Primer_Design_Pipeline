import os
import glob
import subprocess
from Bio import SeqIO


def init(multifasta, directory, config, target, reference):
    combined_seqs = combine_seqs(directory)
    target_db, non_target_db = make_hits_database(combined_seqs, target)
    Constants.assign_values(combined_seqs, config, target, reference, target_db, non_target_db)
    targets = split_multifasta(multifasta)
    return targets


class Constants:

    combined_seqs = None
    config_file = None
    target_list = None
    reference_fasta = None
    target_db = None
    non_target_db = None

    @classmethod
    def assign_values(cls, combined, config, target, reference, target_db, non_target_db):
        cls.combined_seqs = combined
        cls.config_file = config
        cls.target_list = target
        cls.reference_fasta = reference
        cls.target_db = target_db
        cls.non_target_db = non_target_db


class FileNames:

    muscle_input = "tmp_muscle_in.fasta"
    muscle_output = "tmp_muscle_out.fasta"
    conflict_blast_input = "alignment_blast_in.fasta"
    neben_output = "amplicon.txt"
    target_db = "target_database.seqs"
    non_target_db = "non_target_database.seqs"
    target_blast = "target_blast.out"
    non_target_blast = "non_target_blast.out"
    primer3_output = "primer3_out.txt"
    modified_config_file = "config_modified.txt"
    ordering_info = "ordering_info.txt"
    combined_seqs = "combined.seqs"
    
    


def split_multifasta(multifasta):
    out = set()
    with open(multifasta) as infile:
        for record in SeqIO.parse(infile, "fasta"):
            outfile_name = str(record.id) + ".fasta"
            out.add(outfile_name)
            with open(outfile_name, "w") as outfile:
                outfile.write(">{}\n".format(str(record.id)))
                outfile.write(str(record.seq) + "\n")
    return out


def combine_seqs(directory):

    if os.path.isfile("combined.seqs"):
        return os.path.abspath("combined.seqs")

    with open("combined.seqs", "w") as f:

        for file in glob.glob(os.path.join(directory, "*.fasta")):
            f.write(
                ">{}\n".format(os.path.basename(file).replace(".fasta", "")))

            with open(file) as temp_file:
                for record in SeqIO.parse(temp_file, "fasta"):
                    f.write(str(record.seq) + "\n")

    subprocess.run("makeblastdb -in combined.seqs -dbtype nucl > /dev/null 2>&1",
                   shell=True)

    return os.path.abspath("combined.seqs")


def make_hits_database(combined_seqs, target_list):
    targets = set()
    with open(target_list, "rU") as f:
        for line in f:
            targets.add(line.replace("\n", ""))

    with open(FileNames.target_db, "w") as t, open(FileNames.non_target_db, "w") as nt, open(combined_seqs) as f:
        for record in SeqIO.parse(f, "fasta"):
            if str(record.id) in targets:
                t.write(">{}\n{}\n".format(str(record.id), str(record.seq)))
            else:
                nt.write(">{}\n{}\n".format(str(record.id), str(record.seq)))

    subprocess.run("makeblastdb -in {} -dbtype nucl > /dev/null 2>&1".format(FileNames.target_db), shell=True)
    subprocess.run("makeblastdb -in {} -dbtype nucl > /dev/null 2>&1".format(FileNames.non_target_db), shell=True)

    return os.path.abspath(FileNames.target_db), os.path.abspath(FileNames.non_target_db)
