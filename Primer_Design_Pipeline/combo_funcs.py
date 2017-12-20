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
            if forward.value < reverse.value:
                combo_range = abs(forward.value - reverse.value)
                if lower <= combo_range <= upper:
                    temp_combo = Combo(forward, reverse)
                    temp_combo.set_amplicon(upper)
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
        The score is calculated by sorting the list Primers by each attribute.
        After each sort, the position of each Primer object within the list is
        added to the Primer object's score. For example, the primer with the
        fewest degens will have 0 added to its score, while the primer with the
        most degens will have `len(primers)` added to its score.

        If a combo object doesn't have an amplicon, 1000 is added to its score.

    """
    for primer in primers:
        primer.score = 0

    for combo in combos:
        combo.score = 0
    
    def _add_score():
        for i, primer in enumerate(primers):
            primer.score += i

    # num degens
    items = set()
    for primer in primers:
        items.add(primer.num_degens)
    items_list = list(items)
    items_list.sort()
    for primer in primers:
        primer.score += items_list.index(primer.num_degens)


    # Max mis-hit ID
    items = set()
    for primer in primers:
        items.add(primer.max_mis_hit[0])
    items_list = list(items)
    items_list.sort()
    for primer in primers:
        primer.score += items_list.index(primer.max_mis_hit[0])


    # Max mis-hit length
    items = set()
    for primer in primers:
        items.add(primer.max_mis_hit[1])
    items_list = list(items)
    items_list.sort()
    for primer in primers:
        primer.score += items_list.index(primer.max_mis_hit[1])


    # Max non-target hit ID
    items = set()
    for primer in primers:
        items.add(primer.max_non_target_hit[0])
    items_list = list(items)
    items_list.sort()
    for primer in primers:
        primer.score += items_list.index(primer.max_non_target_hit[0])


    # Max non-target hit length
    items = set()
    for primer in primers:
        items.add(primer.max_non_target_hit[1])
    items_list = list(items)
    items_list.sort()
    for primer in primers:
        primer.score += items_list.index(primer.max_non_target_hit[1])


    # Number of amplicons in targets
    items = set()
    for combo in combo:
        items.add(combo.target_amplicons)
    items_list = list(items)
    items_list.sort(reverse=True) # More target amplicons is better
    for combo in combos:
        combo.score += items_list.index(combo.target_amplicons)
    

    # Number of amplicons in non-targets
    items = set()
    for combo in combo:
        items.add(combo.non_target_amplicons)
    items_list = list(items)
    items_list.sort()
    for combo in combos:
        combo.score += items_list.index(combo.non_target_amplicons)


    '''
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
    '''

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
                          "\t# target amplicons,\t# non-target amplicons,"
                          "\tScore\n\n\n")

        outfile.write("\n======================================================"
                      "========================================================"
                      "===============================\n")
        outfile.write(combos[0].target + "\n")
        outfile.write("========================================================"
                      "========================================================"
                      "=============================\n\n")

        format_combos(combos)

        for combo in combos:

            outfile.write("{}.   {}\n".format(combo.id, combo.name))
            outfile.write("----------------------------------------------------"
                          "----------------------------------------------------"
                          "\n")

            outfile.write("{}{}\n".format(combo.forward.display_info, combo.display_info))
            outfile.write("{}{}\n".format(combo.reverse.display_info, combo.display_info))
            outfile.write("\n\n")


def format_combos(combos):

    order = ["name", "mis_hit", "non_target_hit", "num_degens", "tm",
             "amp_len", "target_amp", "non_target_amp", "score"]
    attributes = {}

    for item in order:
        attributes[item] = []

    for combo in combos:
        for primer in [combo.forward, combo.reverse]:
            attributes["name"].append(primer.name)
            attributes["mis_hit"].append(str(primer.max_mis_hit))
            attributes["non_target_hit"].append(str(primer.max_non_target_hit))
            attributes["num_degens"].append(str(primer.num_degens))
            attributes["tm"].append(str(primer.tm))
        attributes["amp_len"].append(str(combo.amp_len))
        attributes["target_amp"].append(str(combo.target_amplicons))
        attributes["non_target_amp"].append(str(combo.non_target_amplicons))
        attributes["score"].append(str(combo.score))

    lengths = []
    for item in order:
        lengths.append(max([len(attribute) for attribute in attributes[item]]) + 4)

    def _format_string(length, attribute):
        out = attribute
        while len(out) < length:
            out += " "
        return out

    for combo in combos:
        combo_attributes = [str(combo.amp_len), str(combo.target_amplicons),
                            str(combo.non_target_amplicons), str(combo.score)]
        for primer in [combo.forward, combo.reverse]:
            attributes = [primer.name, str(primer.max_mis_hit),
                          str(primer.max_non_target_hit),
                          str(primer.num_degens), str(primer.tm)]

            primer_info = ""
            for index, attribute in enumerate(attributes):
                primer_info += _format_string(lengths[index], attribute)

            primer.display_info = primer_info

        combo_info = ""
        combo_info += _format_string(lengths[-4], combo_attributes[-4])
        combo_info += _format_string(lengths[-3], combo_attributes[-3])
        combo_info += _format_string(lengths[-2], combo_attributes[-2])
        combo_info += _format_string(lengths[-1], combo_attributes[-1])

        combo.display_info = combo_info


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
        self.display_info = None

        self.target_amplicons = 0
        self.non_target_amplicons = 0

    @classmethod
    def increment_combo_id(cls):
        cls.combo_id += 1

    def set_amplicon(self, amp_size):
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
            If no amplicon is found, "None found" is assigned instead.

        """
        subprocess.run("{}/Primer_Design_Pipeline/nv2_linux_64 -m {} "
                       "-f {} -r {} -g {} > {}".format(
                           Constants.project_dir,
                           amp_size,
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
            fields = out.split()
            self.amplicon = fields[1]
            self.amp_len = fields[0]

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

