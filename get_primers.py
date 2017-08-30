import subprocess
import os
import glob
from Bio import SeqIO

#TODO
# Return mis_hits and hits
# > Possibly write them to a file
# Use the hits mis_hits and hits when selecting final primers in targetrate_primers.py
# Delete old files


# Driver
# Takes in primer3_template.txt, target.fasta, lower bound, upper bound
# Returns (seqs, target mis-hits, non-target hits)
#         > All dicts with primer name as key, mis-hits and hits have values of lists of tuples with (ID, length)
def get_primers(config_file, target, directory, target_list, lower, upper):
    """Get primers and primer metadata.

    Get primers using primer3. Then find the mis-hits and non-target hits of
    each primer. Mis-hits occur when a primer hits more than one place
    in a single genome. Non-target hits occur when a primer hits an outgroup.

    Args:
        config_file (str): Path to primer3 config file.
        target (str): Path to target genome.
        directory (str): Path to directory with reference .fasta files.
        target_list (str): Path to file containing list of reference genomes.
        lower (int): lower bound of amplicon size.
        upper (int): upper bound of amplicon size.

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
    modify_input_file(config_file, target, lower, upper)
    subprocess.run("primer3_core -output=primer3_out < config_modified.txt", shell=True)
    primers = get_primers_from_primer3("primer3_out")
    
    run_blast(directory, target_list)
    mis_hits, non_target_hits = get_bad_hits(primers)
    
    os.remove("primer3_out")
    os.remove("config_modified.txt")
    
    return (primers, mis_hits, non_target_hits)


# Modify primer3_template.txt to prepare for primer3_core
def modify_input_file(config_file, target, lower, upper):

    # Get marker name and sequence of target.fasta
    with open(target, "U") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            marker_name = record.id
            sequence = str(record.seq)

    # Modify primer3_template.txt and save new version as config_modified.txt
    with open(config_file, "U") as config_in, open("config_modified.txt", "w") as config_out:
        for line in config_in:
            if line.startswith("SEQUENCE_ID"):
                config_out.write("SEQUENCE_ID={}\n".format(marker_name))
            elif line.startswith("SEQUENCE_TEMPLATE"):
                config_out.write("SEQUENCE_TEMPLATE={}\n".format(sequence))
            elif line.startswith("PRIMER_PRODUCT_SIZE_RANGE"):
                config_out.write("PRIMER_PRODUCT_SIZE_RANGE={}-{}\n".format(lower, upper))
            elif line.startswith("PRIMER_NUM_RETURN"):
                pass
            else:
                config_out.write(line)


# Get the primers from primer3_core output
def get_primers_from_primer3(primer3_output):

    # Final output = YY_forward: sequence
    nums = {}   # key = PRIMER_LEFT_X, value = YY
    seqs = {}   # key = PRIMER_LEFT_X, value = sequence

    with open(primer3_output, "U") as infile:
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
    out = set()

    with open("alignment_blast_in.fasta", "w") as f:

        for original_name in seqs:
            primer = Primer(original_name,
                            seqs[original_name],
                            nums[original_name])

            if primer not in out:
                out.add(primer)
                f.write(">{}\n{}\n".format(primer.name, primer.sequence))

    return out


def run_blast(directory, target_list):

    targets = set()
    with open(target_list, "rU") as f:
        for line in f:
            targets.add(line.replace("\n", ""))

    with open("target_database.seqs", "w") as t, open("non_target_database.seqs", "w") as nt, open("combined.seqs") as f:
        for record in SeqIO.parse(f, "fasta"):
            if str(record.id) in targets:
                t.write(">{}\n{}\n".format(str(record.id), str(record.seq)))
            else:
                nt.write(">{}\n{}\n".format(str(record.id), str(record.seq)))

    subprocess.run("makeblastdb -in target_database.seqs -dbtype nucl > /dev/null 2>&1", shell=True)
    subprocess.run("makeblastdb -in non_target_database.seqs -dbtype nucl > /dev/null 2>&1", shell=True)

    subprocess.run("blastn -task blastn -query alignment_blast_in.fasta -db target_database.seqs -num_alignments 2000 -outfmt 6 -out target_blast.out -evalue 10", shell=True)
    subprocess.run("blastn -task blastn -query alignment_blast_in.fasta -db non_target_database.seqs -num_alignments 2000 -outfmt 6 -out non_target_blast.out -evalue 10", shell=True)


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

    with open("target_blast.out", "rU") as f:
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

    with open("non_target_blast.out", "rU") as f:
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
        elif "RIGHT" in original_name:
            self.name = "{}_reverse".format(value)
        
        self.sequence = sequence
        self.value = int(value)
        self.orientation = self.name.split("_")[1]
        self.length = len(sequence)
        self.num_degens = self.count_degens()

        # To be determined at a later time during the program
        self.max_mis_hit = None
        self.max_non_target_hit = None
        self.tm = None, None, None
        self.ordering_info = None

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
        
        lst = [
            target,
            primer_name,
            combined_name,
            final_name,
            sequence_tail,
            to_order,
            self.tm[2],   #Tm
            amplicon,   #amplicon
            len(amplicon),  #amplicon length
            len(amplicon) + 36 # amplicon length + UT length
            ]
        self.ordering_info = ";".join(lst)

    def __hash__(self):
        return hash(self.name + self.sequence)

    def __eq__(self, other):
        if isinstance(other, Primer):
            if self.name == other.name and self.sequence == other.sequence:
                return True
        return False

    def __str__(self):
        return "{}: {}".format(self.name, self.sequence)
