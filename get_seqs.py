import subprocess
import os
import glob
from Bio import SeqIO

#TODO
# Return mis_hits and hits
# > Possibly write them to a file
# Use the hits mis_hits and hits when selecting final primers in generate_primers.py
# Delete old files


# Driver
# Takes in primer3_template.txt, target.fasta, lower bound, upper bound
# Returns (seqs, target mis-hits, non-target hits)
#         > All dicts with primer name as key, mis-hits and hits have values of lists of tuples with (ID, length)
def get_seqs(config_file, gene, directory, target_list, lower, upper):

    modify_input_file(config_file, gene, lower, upper)
    subprocess.run("primer3_core -output=primer3_out < config_modified.txt", shell=True)
    primers = parse_primer3_output("primer3_out")
    
    run_blast(directory, target_list)
    mis_hits, hits = parse_blast_output(primers)
    
    os.remove("primer3_out")
    os.remove("config_modified.txt")
    
    return (primers, mis_hits, hits)




# Modify primer3_template.txt to prepare for primer3_core
def modify_input_file(config_file, gene, lower, upper):

    # Get marker name and sequence of target.fasta
    with open(gene, "U") as infile:
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




# Get the seqs from primer3_core output
def parse_primer3_output(primer3_output):

    # Final output = forward_YY: sequence
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
    out = {}

    with open("alignment_blast_in.fasta", "w") as f:
        for primer in nums:
            
            if "LEFT" in primer:
                #key = "forward_{}".format(nums[primer])
                key = "{}_forward".format(nums[primer])
            else:
                #key = "reverse_{}".format(nums[primer])
                key = "{}_reverse".format(nums[primer])

            if key not in out:
                out[key] = seqs[primer]
                f.write(">{}\n{}\n".format(key, seqs[primer]))

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



def parse_blast_output(primers):
    
    good_hits = {}  # Key = seq name, Value = set of hit genomes
    target_mis_hits = {} # Key = seq name, Value = list of tuples: (mishit ID, mishit length)
    non_target_hits = {}    # Key = seq name, Value = tuple: (hit ID, hit length)


    for primer in primers:
        target_mis_hits[primer] = []
        non_target_hits[primer] = []

    max_mis_hit = (0, 0)
    max_non_target_hit = (0, 0)

    with open("target_blast.out", "rU") as f:
        for line in f:
            fields = line.split()

            # If sequence is in good hits, look for duplicates
            if fields[0] in good_hits:

                # If sequence has hit the same genome in more than one spot                
                if fields[1] in good_hits[fields[0]]:

                    # Append the mis-hit ID and length to target_mis_hits
                    target_mis_hits[fields[0]].append((fields[2], fields[3]))

                    # Keep track of largest mis-hit with ID of at least 90 (subject to change?)
                    if int(fields[3]) > int(max_mis_hit[1]) and float(fields[2]) > 90:
                        max_mis_hit = (float(fields[2]), int(fields[3]))

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
            if int(fields[3]) > int(max_non_target_hit[1]) and float(fields[2]) > 90:
                max_non_target_hit = (float(fields[2]), int(fields[3]))

            
    return ((max_mis_hit, target_mis_hits), (max_non_target_hit, non_target_hits))



if __name__ == "__main__":
    import sys

    args = sys.argv

    print(get_seqs(args[1], args[2], args[3], args[4], args[5], args[6]))

    

    

    
