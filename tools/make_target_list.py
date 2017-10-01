"""The purpose of this file is to provide the user with the ability to create
a target list, which is used to make the combined BLAST database, the target
BLAST database, and the non-target BLAST database.

It should be used after all reference .fasta files are put in a directory, but
before the outgroups are put in. It will create a file called "target_list.txt"
that has the name of each reference genome, one per line.

Pass it in to primer_design_pipeline.py using the -g flag.

"""
import glob
import os
import argparse


def make_target_list(directory):
    """Makes a target_list.

    Makes a file that has the name of each genome in the reference directory,
    one genome name per line.

    Args:
        directory (str): Path to reference directory.

    Returns:
        None.

        Outputs file called "target_list.txt".

    """
    os.chdir(directory)
    with open("target_list.txt", "w") as outfile:
        for file in glob.glob("*.fasta"):
            outfile.write(os.path.splitext(os.path.basename(file))[0] + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="make_target_list.py")

    parser.add_argument("-d", "--directory", help="[REQUIRED] Path to directory"
                        " containing reference fastas", required=True)

    args = parser.parse_args()

    make_target_list(args.directory)
