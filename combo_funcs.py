import subprocess
from setup import *

def get_combos(primers, lower, upper, project_dir):
    out = []

    forwards = []
    reverses = []

    for primer in primers:
        if primer.orientation == "forward":
            forwards.append(primer)
        elif primer.orientation == "reverse":
            reverses.append(primer)

    for f in forwards:
        for r in reverses:
            combo_range = abs(f.value - r.value)
            if lower <= combo_range <= upper:
                temp_combo = Combo(f, r)
                temp_combo.set_amplicon(project_dir)
                out.append(temp_combo)

    return out


def score_combos(primers, combos):
    def _add_score():
        for i in range(len(primers)):
            primers[i].score += i

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

    with open(outfile_name, "a") as outfile:

        outfile.write("Name\tMax mis-hit\tMax non-target hit\t# degens"
                      "\tsequence\t[Min,Max,Avg] Tm\tScore\n\n")
        
        outfile.write("========================================================"
                      "========================================================"
                      "=============================\n")
        outfile.write(combos[0].target + "\n")
        outfile.write("========================================================"
                      "========================================================"
                      "=============================\n\n")

        for combo in combos:

            outfile.write(combo.name + "\n")
            outfile.write("----------------------------------------------------"
                          "----------------------------------------------------"
                          "\n")

            outfile.write("{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\t{}\n".format(
                combo.forward.name,
                combo.forward.max_mis_hit,
                combo.forward.max_non_target_hit,
                combo.forward.num_degens,
                combo.forward.sequence,
                combo.forward.tm,
                combo.score))

            outfile.write("{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\t{}\n".format(
                combo.reverse.name,
                combo.reverse.max_mis_hit,
                combo.reverse.max_non_target_hit,
                combo.reverse.num_degens,
                combo.reverse.sequence,
                combo.reverse.tm,
                combo.score))

            outfile.write("\n\n")


def output_ordering_info(combos):
    with open(FileNames.ordering_info, "w") as outfile:
        for combo in combos:
            outfile.write(combo.name + "\n")
            outfile.write("-------------------------------------------------\n")
            outfile.write(combo.forward.ordering_info + "\n")
            outfile.write(combo.reverse.ordering_info + "\n")
            outfile.write("{}, {}\n".format(combo.target, combo.amplicon))
            outfile.write("\n\n")


def choose_best_combos(combos):
    combos.sort(key=lambda combo: combo.score)
    return combos[0:3]      


class Combo:

    def __init__(self, forward_primer, reverse_primer):
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

    def set_amplicon(self, project_dir):

        subprocess.run("{}/neben_linux_64 -max 500 "
                       "--primers {}:{} {} > {}".format(
                           project_dir,
                           self.forward.sequence,
                           self.reverse.sequence,
                           Constants.reference_fasta,
                           FileNames.neben_output),
                       shell=True)
        
        with open(FileNames.neben_output, "rU") as f:
            out = f.readline()

        if out == "":
            self.amplicon = "None found"
        else:
            self.amplicon = out.split()[3]
