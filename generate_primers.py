# Should ignore_percent be 90% or 10% if you want to ignore 10% of SNPs before considering degens
# Right now it's the former

def generate_primers(sequences, genomes, ignore_percent):

    if ignore_percent > 1:
        ignore_percent /= 100

    # Key = sequence name
    # Val = Primer with amiguity codes
    out = {}

    for name in sequences:

        value, orientation = name.split("_")
        value = int(value)
        length = len(sequences[name])


        if orientation == "reverse":
            sequences[name] = reverse_complement(sequences[name])
        
        degens = [[sequences[name][i]] for i in range(length)]
        
        
        '''
        for genome in genomes:
            
            for i in range(length):
                if orientation == "forward":
                    base = genome[value + i - 1]
                elif orientation == "reverse":
                    base = genome[value - length + i]

                if sequences[name][i] != base:
                    if base != "-":
                        if base not in degens[i]:
                            degens[i].append(base)
          '''

        for i in range(length):

            ignore = 0

            for genome in genomes:
                if orientation == "forward":
                    base = genome[value + i - 1]
                elif orientation == "reverse":
                    base = genome[value - length + i]
                
                if sequences[name][i] != base:
                    if base != "-":
                        if base not in degens[i]:
                            ignore += 1
                            if (ignore / len(genomes)) > (1 - ignore_percent):
                                degens[i].append(base)


        
        primer = ""
        for degen in degens:
            primer += get_code(sorted(degen))

        if "reverse" in name:
            out[name] = reverse_complement(primer)

        else:
            out[name] = primer

    return out



def reverse_complement(sequence):

    out = ""


    codes = {"A":"T", "C":"G", "G":"C", "T":"A", "R":"Y",
             "Y":"R", "S":"S", "W":"W", "K":"M", "M":"K",
             "B":"V", "V":"B", "D":"H", "H":"D", "N":"N"}


    for i in range(len(sequence)-1, -1, -1):
        try:
            out += codes[sequence[i]]
        except KeyError:
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

    '''
    with open("C:\\Users\\Josh\\Desktop\\primer_design_pipeline\\slurm-5821032.out") as f:
        g = f.readlines()[0]


 
    g = g.replace("[", "")
    g = g.replace("]", "")
    g = g.replace("'", "")
    g = g.replace("\n", "")
    
    g = g.split(", ")


    '''
    
    genomes = []
    genome = ""
    with open("C:\\Users\\Josh\\Desktop\\primer_design_pipeline\\all_concatenated_aligned.fasta", "rU") as f:
        for line in f:
            if line.startswith(">"):
                if genome == "":
                    continue
                if "." in genome:
                    print(genome)
                    print(line)
                genomes.append(genome)
                genome = ""
            else:
                genome += line.replace("\n", "")


    


    x = find_degens(s, genomes)
    print(x)
    normal = "ATCG"
    for item in x:
        count = 0
        for name in x[item]:
            if name not in normal:
                count += 1
        print(item + ": " + str(count))
    
