"""Provides functions and class for Combos.

Combos are combinations of Primer objects that fit within the given amplicon
range.

"""
import subprocess
from .setup import Constants, FileNames


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

        if outfile.tell() == 0:
            outfile.write("Key:\n")
            outfile.write("Name,\tMax mis-hit,\tMax non-target hit,\t# degens,"
                          "\tsequence,\t[Min,Max,Avg] Tm,\tAmplicon length,"
                          "\tScore\n\n\n")

        outfile.write("\n======================================================"
                      "========================================================"
                      "===============================\n")
        outfile.write(combos[0].target + "\n")
        outfile.write("========================================================"
                      "========================================================"
                      "=============================\n\n")

        '''
        def _write_primer(combo, primer):
            if combo.amplicon == "None found":
                amp_len = 0
            else:
                amp_len = len(combo.amplicon)
            outfile.write("{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\n".format(
                primer.name,
                primer.max_mis_hit,
                primer.max_non_target_hit,
                primer.num_degens,
                primer.sequence,
                primer.tm,
                amp_len,
                combo.score))
        '''

        for combo in combos:

            outfile.write("{}. {}\n".format(combo.id, combo.name))
            outfile.write("----------------------------------------------------"
                          "----------------------------------------------------"
                          "\n")

            #_write_primer(combo, combo.forward)
            #_write_primer(combo, combo.reverse)
            outfile.write(combo.forward.display_info)
            outfile.write(combo.reverse.display_info)
            outfile.write("\n\n")


def format_combos(combos):

    '''
    name_list = ["Name"]
    mis_hit_list = ["Max mis-hit"]
    non_target_list = ["Max non-target hit"]
    degen_list = ["# degens"]
    sequence_list = ["sequence"]
    tm_list = ["[Min,Max,Avg] Tm"]
    amp_list = ["Amplicon length"]
    score_list = ["Score"]

    for combo in combos:
        name_list.append(combo.forward.name)
        name_list.append(combo.reverse.name)

        mis_hit_list.append(str(combo.forward.mis_hit))
        mis_hit_list.append(str(combo.reverse.mis_hit))

        non_target_list.append(str(combo.forward.non_target_hit))
        non_target_list.append(str(combo.reverse.non_target_hit))

        degen_list.append(str(combo.forward.num_degens))
        degen_list.append(str(combo.reverse.num_degens))

        tm_list.append(str(combo.forward.tm))
        tm_list.append(str(combo.reverse.tm))

        if combo.amplicon == "None found":
            amp_list = 0
        else:
            amp_list = len(combo.amplicon)
        amp_list.append(str(amp_list))

        score_list.append(str(combo.score))

    name_len = _get_max_len(name_list)
    mis_hit_len = _get_max_len(mis_hit_list)
    non_target_len = _get_max_len(non_target_list)
    degen_len = _get_max_len(degen_list)
    sequence_len = _get_max_len(sequence_list)
    tm_len = _get_max_len(tm_list)
    amp_len = _get_max_len(amp_list)
    score_len = _get_max_len(score_list)
    '''

    attribute_list = []
    for combo in combos:
        attribute_list.append(combo.forward.name)
        attribute_list.append(combo.reverse.name)
        attribute_list.append(str(combo.forward.max_mis_hit))
        attribute_list.append(str(combo.reverse.max_mis_hit))
        attribute_list.append(str(combo.forward.max_non_target_hit))
        attribute_list.append(str(combo.reverse.max_non_target_hit))
        attribute_list.append(str(combo.forward.num_degens))
        attribute_list.append(str(combo.reverse.num_degens))
        attribute_list.append(str(combo.forward.tm))
        attribute_list.append(str(combo.reverse.tm))
        attribute_list.append(str(combo.amp_len))
        attribute_list.append(str(combo.score))

    max_len = max(len(attribute) for attribute in attribute_list) + 2

    def _format_string(attribute):
        out = attribute
        while len(out) < max_len:
            out += " "
        return out

    for combo in combos:
        for primer in [combo.forward, combo.reverse]:
            attributes = [primer.name, primer.max_mis_hit,
                          primer.max_non_target_hit, primer.num_degens,
                          str(primer.tm), str(combo.amp_len),
                          str(combo.score)]
            out = ""
            for attribute in attributes:
                out += _format_string(attribute)
            primer.display_info = out

    #TODO
    # Add a display_info instance variable to Primer class.
    #   > Just one string that is formattted for that primer
    # In output_primers, just write forward.display_info \n reverse.display_info


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

    combo_id = 1

    def __init__(self, forward_primer, reverse_primer):
        """Creates a Combo object sans amplicon, ordering info, and score.

        Args:
            forward_primer (:obj:`Primer`): Primer object with forward
                orientation.
            reverse_primer (:obj:`Primer`): Primer object with reverse
                orientation.

        """
        self.id = Combo.combo_id
        Combo.increment_combo_id()
        
        self.forward = forward_primer
        self.reverse = reverse_primer
        self.name = "{} - {}".format(self.forward.name, self.reverse.name)
        self.amplicon = None
        self.amp_len = None
        self.target = None
        self.primer_name = None
        self.combined_name = None
        self.sequence = None
        self.target = None
        self.score = 0

    @classmethod
    def increment_combo_id(cls):
        cls.combo_id += 1

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
        subprocess.run("{}/Primer_Design_Pipeline/neben_linux_64 -max 500 "
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
            self.amp_len = 0
        else:
            self.amplicon = out.split()[3]
            self.amp_len = len(self.amplicon)

    def get_order_info(self):
        """Gets the ordering info for the forward and reverse primers.

        Gets the ordering info for the forward and reverse primers in the form
        of a list. Can be used to write the ordering info to a csv file.

        Args:
            None.

        Returns:
            :obj:`list` of :obj:`list` of :obj:`str`: A list of the ordering
            information for the forward and reverse primers in the combo.

        """
        out = []
        for primer in [self.forward, self.reverse]:
            out.append(["", self.target, primer.ordering_info["primer_name"],
                        primer.ordering_info["combined_name"], primer.sequence,
                        primer.ordering_info["final_name"],
                        primer.ordering_info["sequence_tail"],
                        primer.ordering_info["to_order"],
                        primer.ordering_info["tm"],
                        self.amplicon, str(len(self.amplicon)),
                        str(len(self.amplicon)+36), "", "", "", ""])
        return out

    def get_amplicon_info(self):
        return [self.target, self.amplicon]
