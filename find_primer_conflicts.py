import os
import subprocess

# Driver
# Takes in fasta file with all primers
def find_primer_conflicts(primer_fasta):

    primer_path = os.path.abspath(primer_fasta)
    subprocess.run("makeblastdb -in {} -dbtype nucl > /dev/null 2>&1".format(primer_path), shell=True)
    subprocess.run("blastn -task blastn -query {} -db {} -perc_identity 100 -dust no -evalue 20 -word_size 4 -outfmt 6 -out primer_conflicts_blast.out > /dev/null 2>&1".format(primer_path, primer_path), shell=True)

    out = parse_blast_output("primer_conflicts_blast.out")
    os.remove("primer_conflicts_blast.out")

    output_conflicts(out)

    return out



def parse_blast_output(output):

    # Key = alignment length
    # Value = list of conflicts that occur at that alignment length
    #         Each conflict is a list of length 2 where item 0 conflicts with item 1
    conflicts = {}

    # Lower keys should contain the values of higher keys, but don't.
    # This needs to be accounted for in whatever script calls this one.
    # E.g. conflicts[9] should contain both the values of conflicts[9] and of conflicts[10].
    # conflicts[8] should contain conflicts[8], conflicts[9], and conflicts[10].  Etc. etc.
    

    # If there are no conflicts at a certain length, conflicts[length] == []
    for i in range(1, 11):
        conflicts[i] = []

    with open(output, "rU") as f:

        for line in f:
            fields = line.strip().split()

            # What is fields[5]?
            # fields[3] = length, fields[0] = primer #1, fields[1] = primer #2
            if int(fields[3]) <= 10 and int(fields[5]) == 0:
                conflicts[int(fields[3])].append([fields[0], fields[1]])

    print(conflicts)
    return conflicts



def output_conflicts(conflicts):

    '''
    with open("primer_conflicts.txt", "w") as out:
        for length in conflicts:
            if len(conflicts[length]) > 0:
                for i in range(10, length-1, -1):
                    out.write("Problematic primers with an alignment length of {}\n".format(length))
                    out.write("{}\t{}\n\n".format(conflicts[i][0], conflicts[i][1]))
    '''
    pass


def blast_all_primers(primer_fasta, combined_seqs):
    subprocess.run("blastall -p blastn -i {} -d {} -m 8 -e 10 -o primers.blast.out".format(primer_fasta, combined_seqs), shell=True)
