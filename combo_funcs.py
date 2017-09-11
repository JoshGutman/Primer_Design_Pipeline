"""Provides functions and class for Combos.

Combos are combinations of Primer objects that fit within the given amplicon
range.

"""
import subprocess
from setup import Constants, FileNames


def get_combos(primers, lower, upper):
    """Gets possible combos within given amplicon range.

    Args:
        primers (:obj:`list` of :obj:`Primer`): List of Primer objects.
        lower (int): Lower bound of amplicon size.
        upper (int): Upper bound of amplicon size.

    Returns:
        :obj:`list` of :obj:`Combo`: List of Combo objects.

    """
    out = []
    forwards = []
    reverses = []

    for primer in primers:
        if primer.orientation == "forward":
            forwards.append(primer)
        elif primer.orientation == "reverse":
            reverses.append(primer)

    # Find combinations of Primers that fall within lower - upper
    for forward in forwards:
        for reverse in reverses:
            combo_range = abs(forward.value - reverse.value)
            if lower <= combo_range <= upper:
                temp_combo = Combo(forward, reverse)
                temp_combo.set_amplicon()
                out.append(temp_combo)
    return out


def score_combos(primers, combos):
    """Set Combo scores.

    Set Combo scores based on the number of degens and the max mis-hit and
    non-target hit IDs and lengths.

    Args:
        primers (:obj:`list` of :obj:`Primer`): List of Primer objects.
        combos (:obj:`list` of :obj:`Combo`): List of Combo objects.

    Returns:
        None.

        Mutates the Primer.score and Combo.score instance variables.

    Notes:
        The score is calculated by sorting the list Primers by each factor.
        After each sort, the position of each Primer object within the list is
        added to the Primer object's score. For example, the primer with the
        fewest degens will have 0 added to its score, while the primer with the
        most degens will have `len(primers)` added to its score.

        If a combo object doesn't have an amplicon, 1000 is added to its score.

    """
    def _add_score():
        for i, primer in enumerate(primers):
            primer.score += i

    # Number of degens
    primers.sort(key=lambda primer: primer.num_degens)
    _add_score()

    # Max mis-hit ID
    primers.sort(key=lambda primer: primer.max_mis_hit[0])
    _add_score()

    # Max mis-hit length
    primers.sort(key=lambda primer: primer.max_mis_hit[1])
    _add_score()

    # Max non-target hit ID
    primers.sort(key=lambda primer: primer.max_non_target_hit[0])
    _add_score()

    # Max non-target hit length
    primers.sort(key=lambda primer: primer.max_non_target_hit[1])
    _add_score()

    for combo in combos:
        combo.score += combo.forward.score + combo.reverse.score
        if combo.amplicon == "None found":
            combo.score += 1000


def output_combos(combos, outfile_name):
    """Outputs Combo information for comparison.

    Args:
        combos (:obj:`list` of :obj:`Combo`): List of Combo objects.
        outfile_name (str): Name of file to be written to.

    Returns:
        None.

        Outputs a file defined by outfile_name.

    """
    with open(outfile_name, "a") as outfile:

        outfile.write("Name\tMax mis-hit\tMax non-target hit\t# degens"
                      "\tsequence\t[Min,Max,Avg] Tm\tAmplicon length\tScore\n\n")

        outfile.write("========================================================"
                      "========================================================"
                      "=============================\n")
        outfile.write(combos[0].target + "\n")
        outfile.write("========================================================"
                      "========================================================"
                      "=============================\n\n")

        def _write_primer(combo, primer):
            outfile.write("{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\t{}\n".format(
                primer.name,
                primer.max_mis_hit,
                primer.max_non_target_hit,
                primer.num_degens,
                primer.sequence,
                primer.tm,
                len(combo.amplicon),
                combo.score))

        for combo in combos:

            outfile.write(combo.name + "\n")
            outfile.write("----------------------------------------------------"
                          "----------------------------------------------------"
                          "\n")

            _write_primer(combo, combo.forward)
            _write_primer(combo, combo.reverse)
            outfile.write("\n\n")


def output_ordering_info(combos):
    """Writes ordering info of combos to a file.

    Args:
        combos (:obj:`list` of :obj:`Combo`): List of Combo objects.

    Returns:
        None.

        Outputs file defined in FileNames.ordering_info.

    """
    with open(FileNames.ordering_info, "w") as outfile:
        for combo in combos:
            outfile.write(combo.name + "\n")
            outfile.write("-------------------------------------------------\n")
            outfile.write(combo.forward.ordering_info + "\n\n")
            outfile.write(combo.reverse.ordering_info + "\n\n")
            outfile.write("{}, {}\n".format(combo.target, combo.amplicon))
            outfile.write("\n\n\n")


def choose_best_combos(combos):
    """Choose the three combos with highest scores.

    Args:
        combos (:obj:`list` of :obj:`Combo`): List of Combo objects.

    Returns:
        :obj:`list` of :obj:`Combo`: List of 3 Combos with highest scores.

    """
    combos.sort(key=lambda combo: combo.score)
    return combos[0:3]


class Combo:
    """Combos consist of two Primer objects and related attributes.

    Realted attributes include name, amplicon, target genome, ordering
    information, and score.

    """

    def __init__(self, forward_primer, reverse_primer):
        """Creates a Combo object sans amplicon, ordering info, and score.

        Args:
            forward_primer (:obj:`Primer`): Primer object with forward
                orientation.
            reverse_primer (:obj:`Primer`): Primer object with reverse
                orientation.

        """
        self.forward = forward_primer
        self.reverse = reverse_primer
        self.name = "{} - {}".format(self.forward.name, self.reverse.name)
        self.amplicon = None
        self.target = None
        self.primer_name = None
        self.combined_name = None
        self.sequence = None
        self.target = None
        self.score = 0

    def set_amplicon(self):
        """Sets the amplicon of a Combo object.

        Runs neben to find the amplicon(s) of a Combo object. The output of
        neben is redirected to a file defined in FileNames.neben_output. The
        output is then read and assigned to Combo.amplicon.

        Args:
            None.

        Returns:
            None.

            Mutates the Combo.amplicon instance variable.

        Notes:
            The max size of an amplicon is set to 500.

            If no amplicon is found, "None found" is assigned instead.

        """
        subprocess.run("{}/neben_linux_64 -max 500 "
                       "--primers {}:{} {} > {}".format(
                           Constants.project_dir,
                           self.forward.sequence,
                           self.reverse.sequence,
                           Constants.reference_fasta,
                           FileNames.neben_output),
                       shell=True)

        with open(FileNames.neben_output, "rU") as infile:
            out = infile.readline()

        if out == "":
            self.amplicon = "None found"
        else:
            self.amplicon = out.split()[3]
