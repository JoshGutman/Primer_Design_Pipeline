#TODO
# Something wrong with reverse primer generation
#   > base = genome[value - length + i - 1] incorrect?
#   > Another reverse_complement needed somewhere?

def find_degens(sequences, genomes):

    # Key = sequence name
    # Val = Primer with amiguity codes
    out = {}

    for seq in sequences:

        orientation, value = seq.split("_")
        value = int(value)
        sequence = sequences[seq]
        length = len(sequence)
        
        degens = [[sequence[i]] for i in range(length)]
        
        if orientation == "reverse":
            sequences[seq] = reverse_complement(sequences[seq])


        for genome in genomes:
            
            for i in range(length):
                if orientation == "forward":
                    base = genome[value + i - 1]
                elif orientation == "reverse":
                    base = genome[value - length + i - 1]
                    
                if sequence[i] != base:
                    if base != "-":
                        if base not in degens[i]:
                            degens[i].append(base)

        primer = ""
        for degen in degens:
            primer += get_code(sorted(degen))

        out[seq] = primer

    return out



def reverse_complement(sequence):

    out = ""

    for i in range(len(sequence) - 1, -1, -1):

        if sequence[i] == 'A':
            out += 'T'
        elif sequence[i] == 'C':
            out += 'G'
        elif sequence[i] == 'G':
            out += 'C'
        elif sequence[i] == 'T':
            out += 'A'
        else:
            raise ValueError("Non-ACGT char encountered when attempting to compute reverse complement\nSequence: {}".format(sequence))

    return out
                  




def get_code(bases):

    if len(bases) == 1:
        return "".join(bases)

    elif bases == ["A", "G"]:
      return "R"
    elif bases == ["C", "T"]:
      return "Y"
    elif bases == ["C", "G"]:
      return "S"
    elif bases == ["A", "T"]:
      return "W"
    elif bases == ["G", "T"]:
      return "K"
    elif bases == ["A", "C"]:
      return "M"
    elif bases == ["C", "G", "T"]:
      return "B"
    elif bases == ["A", "G", "T"]:
      return "D"
    elif bases == ["A", "C", "T"]:
      return "H"
    elif bases == ["A", "C", "G"]:
      return "V"
    elif bases == ["A", "C", "G", "T"]:
      return "N"
    else:
        raise ValueError("Unrecognized sequence of bases encountered when attempting to compute the ambiguity code\nSequence: {}".format(bases))

if __name__ == "__main__":

    s = {'forward_167': 'CATCCGCAATGAGAGTAAAT', 'forward_99': 'ATAGTATCAACGCCGACATC', 'reverse_380': 'GGTCTGGATGGTTGAACTTA', 'reverse_244': 'GGCAAAAGCTATTTTCTCAA', 'forward_193': 'TTATCTCTTTCAGGGTCGTG', 'reverse_394': 'ATATCATAACGGGTGGTCTG', 'reverse_291': 'TTTTTCGTGAACTGGAAACT', 'forward_15': 'ATGGTTGGTGTCAGATCTTC'}
    with open("..\\slurm-5821032.out", "rU") as f:
        g = f.readlines()[0]

    
    g = g.replace("[", "")
    g = g.replace("]", "")
    g = g.replace("'", "")
    g = g.replace("\n", "")
    
    g = g.split(", ")

    print(reverse_complement("GGCAAAAGTTATTTCCTCAA"))
    print(find_degens(s, g))
