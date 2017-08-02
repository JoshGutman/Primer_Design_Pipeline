import subprocess
import os
from Bio import SeqIO

# Driver
# Takes in primer3_template.txt, target.fasta, lower bound, upper bound
def get_seqs(config_file, gene, lower, upper):

    modify_input_file(config_file, gene, lower, upper)
    subprocess.run("primer3_core -output=primer3_out < config_modified.txt", shell=True)
    out = parse_primer3_output("primer3_out")
    os.remove("primer3_out")
    os.remove("config_modified.txt")
    return out




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
    for primer in nums:
        
        if "LEFT" in primer:
            key = "forward_{}".format(nums[primer])
        else:
            key = "reverse_{}".format(nums[primer])
            
        out[key] = seqs[primer]

    return out





if __name__ == "__main__":
    import sys

    args = sys.argv

    print(get_seqs(args[1], args[2], args[3], args[4]))

    

    

    
