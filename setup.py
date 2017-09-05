import os
import glob
import subprocess
from Bio import SeqIO


def init(multifasta, directory, config, target, reference, keep):
    combined_seqs = combine_seqs(directory)
    Constants.assign_values(combined_seqs, config, target, reference, keep)
    targets = split_multifasta(multifasta)
    return targets


class Constants:

    combined_seqs = None
    config_file = None
    target_list = None
    reference_fasta = None
    keep_files = None

    @classmethod
    def assign_values(cls, combined, config, target, reference, keep):
        cls.combined_seqs = combined
        cls.config_file = config
        cls.target_list = target
        cls.reference_fasta = reference
        cls.keep_files = keep


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
