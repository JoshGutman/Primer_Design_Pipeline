"""Calculates the degeneracies of a list of Primer objects.

Mutates the passed-in Primer objects.

"""
def get_degens(primers, genomes, ignore_percent):
    """Calculates and sets the degens for a list of primers.

    Given a list of Primer objects, a list of the genomes, and an ignore
    percentage, a new sequence will be created with the degeneracies calculated.

    Args:
        primers (:obj:`list` of :obj:`Primer`): List of Primer objects.
        genomes (:obj:`list` of :obj:`str`): List of genomes.
        ignore_percent(float): A SNP will have to have to occur in at least
            (1 - ignore_percent)% genomes for it to be considered a degen.

    Returns:
        None.

        Mutates the `sequence` attribute of the passed in `Primer`s.

    Notes:
        Driver function.

        When looking for degens, an unique SNP (meaning the base and position)
        will have to occur in at least (1 - ignore_percent)% of genomes before
        it's considered a degen.

    """
    if ignore_percent > 1:
        ignore_percent /= 100

    for primer in primers:

        if primer.orientation == "reverse":
            primer.sequence = _reverse_complement(primer.sequence)

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

                if (primer.sequence[i] != base
                        and base != "-" and
                        base not in degens[i]):

                    snps[base] += 1

                    # Only consider degens if the base occurs in more
                    # than (100 - ignore_percent) of genomes
                    current_percentage = snps[base] / len(genomes)
                    if (current_percentage) > (1 - ignore_percent):
                        degens[i].append(base)

        new_seq = ""
        for degen in degens:
            new_seq += _get_code(sorted(degen))

        if primer.orientation == "reverse":
            new_seq = _reverse_complement(new_seq)

        primer.set_sequence(new_seq)


def _reverse_complement(sequence):
    """Gets the reverse complement of a sequence.

    Args:
        sequence (str): Primer's sequence containing ACGT chars.

    Returns:
        str: The reverse complement of `sequence`.

    Raises:
        ValueError: If sequence contains a non-ACGT char.

    """
    out = ""

    codes = {"A":"T", "C":"G", "G":"C", "T":"A", "R":"Y",
             "Y":"R", "S":"S", "W":"W", "K":"M", "M":"K",
             "B":"V", "V":"B", "D":"H", "H":"D", "N":"N"}

    sequence = sequence.upper()

    for i in range(len(sequence)-1, -1, -1):
        try:
            out += codes[sequence[i]]
        except KeyError:
            raise ValueError("Non-ACGT char encountered when attempting to "
                             "compute reverse complement.\nSequence:"
                             " {}".format(sequence))

    return out


def _get_code(bases):
    """Gets the nucleotide ambiguity code.

    Args:
        bases (:obj:`list` of :obj:`str`): List of bases that occured in a
            position.

    Returns:
        str: Nucleotide ambiguity code.

    Raises:
        ValueError: If a non-ACGT char is found in `bases` or
            if `bases` is not sorted.

    """
    bases_str = "".join(bases)

    degen_dict = {"AG": "R", "CT": "Y", "CG": "S",
                  "AT": "W", "GT": "K", "AC": "M",
                  "CGT": "B", "AGT": "D", "ACT": "H",
                  "ACG": "V", "ACGT": "N"}

    if len(bases) == 1:
        return bases_str

    if bases_str not in degen_dict:
        raise ValueError("Unrecognized sequence of bases encountered when "
                         "attempting to compute the ambiguity code\nSequence:"
                         " {}".format(bases))
    else:
        return degen_dict[bases_str]
