# Should ignore_percent be 90% or 10% if you want to ignore 10% of SNPs before considering degens
# Right now it's the former
import sys

def get_degens(primers, genomes, ignore_percent):

    if ignore_percent > 1:
        ignore_percent /= 100

    for primer in primers:

        if primer.orientation == "reverse":
            primer.sequence = reverse_complement(primer.sequence)

        degens = [[primer.sequence[i]] for i in range(primer.length)]

        for i in range(primer.length):

            snps = {"A": 0,
                    "C": 0,
                    "G": 0,
                    "T": 0}

            for genome in genomes:
                if primer.orientation == "forward":
                    base = genome[primer.value + i - 1]
                elif primer.orientation == "reverse":
                    base = genome[primer.value - primer.length + i]
                else:
                    print("ERROR 0")
                    sys.exit()
                
                if primer.sequence[i] != base:
                    if base != "-":
                        if base not in degens[i]:
                            snps[base] += 1

                            # Only consider degens if the base occurs in more than (100 - ignore_percent) of genomes
                            if (snps[base] / len(genomes)) > (1 - ignore_percent):
                                degens[i].append(base)
        
        new_seq = ""
        for degen in degens:
            new_seq += get_code(sorted(degen))

        if primer.orientation == "reverse":
            new_seq = reverse_complement(new_seq)

        primer.set_sequence(new_seq)


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


