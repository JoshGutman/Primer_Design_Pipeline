import subprocess
import os
import glob
from Bio import SeqIO

from .setup import Constants, FileNames

# Driver
# Takes in primer3_template.txt, target.fasta, lower bound, upper bound
# Returns (seqs, target mis-hits, non-target hits)
#         > All dicts with primer name as key, mis-hits and hits have values of lists of tuples with (ID, length)
def get_primers(config_file, target, directory, lower, upper, temp_args):
    """Get primers and primer metadata.

    Get primers using primer3. Then find the mis-hits and non-target hits of
    each primer. Mis-hits occur when a primer hits more than one place
    in a single genome. Non-target hits occur when a primer hits an outgroup.

    Args:
        config_file (str): Path to primer3 config file.
        target (str): Path to target genome.
        directory (str): Path to directory with reference .fasta files.
        lower (int): Lower bound of amplicon size.
        upper (int): Upper bound of amplicon size.
        temp_args (:obj:`list` of :obj:`float`): List of ol_conc,na_conc,mg_conc

    Returns:
        tuple: Tuple consisting of primers, mis-hits, and non-target hits.

            primers (dict): Keys are the primer name (e.g. 244_reverse), values
            are the primer sequence.

            mis-hits (dict): Keys are the primer name, values are lists of
            tuples, each tuple being a mis-hit ID and length.

            non-target hits (dict): Keys are the primer name, values are lists
            of tuples, each tuple being a non-target hit ID and length.

    Notes:


    """
    modify_input_file(config_file, target, lower, upper, temp_args)
    subprocess.run("primer3_core -output={} < {}".format(
        FileNames.primer3_output, FileNames.modified_config_file),
                   shell=True)
    primers = get_primers_from_primer3()

    run_blast(directory)
    mis_hits, non_target_hits = get_bad_hits(primers)

    return (primers, mis_hits, non_target_hits)


# Modify primer3_template.txt to prepare for primer3_core
def modify_input_file(config_file, target, lower, upper, temp_args):

    # Get marker name and sequence of target.fasta
    with open(target, "U") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            marker_name = record.id
            sequence = str(record.seq)

    # Modify primer3_template.txt and save new version as config_modified.txt
    with open(config_file, "U") as config_in, open(FileNames.modified_config_file, "w") as config_out:
        for line in config_in:
            if line.startswith("SEQUENCE_ID"):
                config_out.write("SEQUENCE_ID={}\n".format(marker_name))
            elif line.startswith("SEQUENCE_TEMPLATE"):
                config_out.write("SEQUENCE_TEMPLATE={}\n".format(sequence))
            elif line.startswith("PRIMER_PRODUCT_SIZE_RANGE"):
                config_out.write("PRIMER_PRODUCT_SIZE_RANGE={}-{}\n".format(lower, upper))

##            elif line.startswith("PRIMER_SALT_MONOVALENT"):
##                config_out.write("PRIMER_SALT_MONOVALENT={}\n".format(temp_args[1]))
##            elif line.startswith("PRIMER_INTERNAL_SALT_MONOVALENT"):
##                config_out.write("PRIMER_INTERNAL_SALT_MONOVALENT={}\n".format(temp_args[1]))
##            elif line.startswith("PRIMER_SALT_DIVALENT"):
##                config_out.write("PRIMER_SALT_DIVALENT={}\n".format(temp_args[2]))
##            elif line.startswith("PRIMER_INTERNAL_SALT_DIVALENT"):
##                config_out.write("PRIMER_INTERNAL_SALT_DIVALENT={}\n".format(temp_args[2]))
##            elif line.startswith("PRIMER_DNA_CONC"):
##                config_out.write("PRIMER_DNA_CONC={}\n".format(int(temp_args[0])*1000))
##            elif line.startswith("PRIMER_INTERNAL_DNA_CONC"):
##                config_out.write("PRIMER_INTERNAL_DNA_CONC={}\n".format(int(temp_args[0])*1000))
##            elif line.startswith("PRIMER_SALT_CORRECTIONS"):
##                config_out.write("PRIMER_SALT_CORRECTIONS=2\n")
##            elif line.startswith("PRIMER_TM_FORMULA"):
##                config_out.write("PRIMER_TM_FORMULA=1\n")
                
            elif line.startswith("PRIMER_NUM_RETURN"):
                pass
            else:
                config_out.write(line)


# Get the primers from primer3_core output
def get_primers_from_primer3():

    # Final output = YY_forward: sequence
    nums = {}   # key = PRIMER_LEFT_X, value = YY
    seqs = {}   # key = PRIMER_LEFT_X, value = sequence

    with open(FileNames.primer3_output, "U") as infile:
        for line in infile:

            # PRIMER_LEFT_X or PRIMER_RIGHT_X
            if "," in line and " " not in line:
                fields = line.split("=")
                nums[fields[0]] = fields[1].split(",")[0]

            # PRIMER_LEFT_X_SEQUENCE or PRIMER_RIGHT_X_SEQUENCE
            elif "_SEQUENCE" in line:
                fields = line.split("=")
                seqs[fields[0].replace("_SEQUENCE", "")] = fields[1].replace("\n", "")

    # Get the primer value and sequence from nums and seqs
    out = []

    with open(FileNames.conflict_blast_input, "w") as f:

        for original_name in seqs:
            primer = Primer(original_name,
                            seqs[original_name],
                            nums[original_name])

            if primer not in out:
                out.append(primer)
                f.write(">{}\n{}\n".format(primer.name, primer.sequence))

    return out


def run_blast(directory):

    subprocess.run("blastn -task blastn -query {} -db {} -num_alignments 2000 -outfmt 6 -out {} -evalue 10".format(FileNames.conflict_blast_input, Constants.target_db, FileNames.target_blast), shell=True)
    subprocess.run("blastn -task blastn -query {} -db {} -num_alignments 2000 -outfmt 6 -out {} -evalue 10".format(FileNames.conflict_blast_input, Constants.non_target_db, FileNames.non_target_blast), shell=True)


def get_bad_hits(primers):

    good_hits = {}  # Key = seq name, Value = set of hit genomes
    mis_hits = {} # Key = seq name, Value = list of tuples: (mishit ID, mishit length)
    non_target_hits = {} # Key = seq name, Value = tuple: (hit ID, hit length)
    max_mis_hits = {}
    max_non_target_hits = {}

    for primer in primers:
        mis_hits[primer.name] = []
        non_target_hits[primer.name] = []

        max_mis_hits[primer.name] = (0, 0)
        max_non_target_hits[primer.name] = (0, 0)

    with open(FileNames.target_blast, "rU") as f:
        for line in f:
            fields = line.split()

            # If sequence is in good hits, look for duplicates
            if fields[0] in good_hits:

                # If sequence has hit the same genome in more than one spot
                if fields[1] in good_hits[fields[0]]:

                    # Append the mis-hit ID and length to mis_hits
                    mis_hits[fields[0]].append((fields[2], fields[3]))

                    # Keep track of largest mis-hit with ID of at least 90 (subject to change?)
                    if int(fields[3]) > int(max_mis_hits[fields[0]][1]) and float(fields[2]) > 90:
                        max_mis_hits[fields[0]] = (float(fields[2]), int(fields[3]))

                # If sequence has hit a genome for the first time, add it to good_hits
                else:
                    good_hits[fields[0]].add(fields[1])

            # If sequence not is not in good_hits, set the initial value to an empty set
            else:
                good_hits[fields[0]] = set()

    with open(FileNames.non_target_blast, "rU") as f:
        for line in f:
            fields = line.split()

            non_target_hits[fields[0]].append((fields[2], fields[3]))

            # Keep track of largest non-target hit with ID of at least 90 (subject to change?)
            if int(fields[3]) > int(max_non_target_hits[fields[0]][1]) and float(fields[2]) > 90:
                max_non_target_hits[fields[0]] = (float(fields[2]), int(fields[3]))

    for primer in primers:
        primer.max_mis_hit = max_mis_hits[primer.name]
        primer.max_non_target_hit = max_non_target_hits[primer.name]

    return mis_hits, non_target_hits


class Primer:

    def __init__(self, original_name, sequence, value):

        if "LEFT" in original_name:
            self.name = "{}_forward".format(value)
            self.orientation = "forward"
        elif "RIGHT" in original_name:
            self.name = "{}_reverse".format(value)
            self.orientation = "reverse"

        self.sequence = sequence
        self.value = int(value)
        self.length = len(sequence)
        self.num_degens = self.count_degens()

        # To be determined at a later time during the program
        self.max_mis_hit = None, None
        self.max_non_target_hit = None, None
        self.tm = None, None, None
        self.ordering_info = None
        self.display_info = None
        self.score = 0

    def set_sequence(self, new_seq):
        self.sequence = new_seq
        self.length = len(self.sequence)
        self.num_degens = self.count_degens()

    def count_degens(self):
        out = 0
        for base in self.sequence:
            if base not in "ACGT":
                out += 1
        return out

    def set_ordering_info(self, target, amplicon):

        self.ordering_info = {}

        ut1 = "ACCCAACTGAATGGAGC"
        ut2 = "ACGCACTTGACTTGTCTTC"

        tails = {"forward": ("UT1", ut1),
                 "reverse": ("UT2", ut2)}

        target = os.path.splitext(target)[0]
        primer_name = str(self.value) + self.orientation[0].upper()
        combined_name = "{}_{}".format(primer_name, target)
        final_name = "{}_{}".format(combined_name, tails[self.orientation][0])
        sequence_tail = tails[self.orientation][1] + self.sequence
        to_order = "{},{}".format(final_name, sequence_tail)

        self.ordering_info["target"] = target
        self.ordering_info["primer_name"] = primer_name
        self.ordering_info["combined_name"] = combined_name
        self.ordering_info["final_name"] = final_name
        self.ordering_info["sequence_tail"] = sequence_tail
        self.ordering_info["to_order"] = to_order
        self.ordering_info["tm"] = str(self.tm[2])

    def __hash__(self):
        return hash(self.name + self.sequence)

    def __eq__(self, other):
        if isinstance(other, Primer):
            if self.name == other.name and self.sequence == other.sequence:
                return True
        return False

    def __str__(self):
        return "{}: {}".format(self.name, self.sequence)
