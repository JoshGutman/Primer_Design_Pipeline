import subprocess
import os
from Bio import SeqIO

def create_seqs(config_file, gene, upper, lower):

    out = []

    # Not needed if only one genome in target.fasta
    with open(os.path.abspath(gene), "U") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            gene_name = record.id
            modify_input_file(config_file, gene, gene_name, upper, lower)
            subprocess.run("primer3_core -output=primer3_out < config_modified.txt", shell=True)
            out.append(parse_primer3_output("primer3_out"))
    return out



def modify_input_file(config_file, gene, gene_name, upper, lower):

    sequence = ""

    # Unneeded if only one genome in target.fasta
    with open(gene, "U") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id == gene_name:
                marker_name = record.id
                sequence += str(record.seq)

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



def parse_primer3_output(primer3_output):

    nums = {}
    seqs = {}

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

    print(create_seqs(args[1], args[2], args[3], args[4]))

    

    

    
