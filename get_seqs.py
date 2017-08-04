import subprocess
import os
import glob
from Bio import SeqIO

'''
ecoli-TARGET_LIST.TXT

'''

# Driver
# Takes in primer3_template.txt, target.fasta, lower bound, upper bound
def get_seqs(config_file, gene, directory, target_list, lower, upper):

    modify_input_file(config_file, gene, lower, upper)
    subprocess.run("primer3_core -output=primer3_out < config_modified.txt", shell=True)
    out = parse_primer3_output("primer3_out")
    get_all_alignments(directory, target_list)
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




def get_all_alignments(directory, target_list):

    targets = set()
    with open(target_list, "rU") as f:
        for line in f:
            targets.add(line.replace("\n", ""))


    with open("target_database.seqs", "w") as t, open("non_target_database.seqs", "w") as nt:
        for file in glob.glob(os.path.join(directory, "*.fasta")):
            name = os.path.basename(file).replace(".fasta", "")

            if name in targets:
                t.write(">{}\n".format(name))
                with open(file) as f:
                    for record in SeqIO.parse(f, "fasta"):
                        t.write(str(record.seq) + "\n")
            else:
                nt.write(">{}\n".format(name))
                with open(file) as f:
                    for record in SeqIO.parse(f, "fasta"):
                        nt.write(str(record.seq) + "\n")
            
    subprocess.run("makeblastdb -in target_database.seqs -dbtype nucl > /dev/null 2>&1", shell=True)
    subprocess.run("makeblastdb -in non_target_database.seqs -dbtype nucl > /dev/null 2>&1", shell=True)

    subprocess.run("blastn -task blastn -query alignment_blast_in.fasta -db target_database.seqs -num_alignments 2000 -outfmt 6 -out target_blast.out -evalue 10", shell=True)
    subprocess.run("blastn -task blastn -query alignment_blast_in.fasta -db non_target_database.seqs -num_alignments 2000 -outfmt 6 -out non_target_blast.out -evalue 10", shell=True)




if __name__ == "__main__":
    import sys

    args = sys.argv

    print(get_seqs(args[1], args[2], args[3], args[4], args[5], args[6]))

    

    

    
