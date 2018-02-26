import subprocess
import argparse
import pickle
import glob
import sys
import os


def database_amplicons(directory, execute, primer_id, forward, reverse, amp_size, target, keep):

    primer = _interpret_primer(primer_id, forward, reverse)
    neben_path = _get_neben_path()

    files = _get_files(directory)
    num_files = len(files)
    groups, species_names = _split_files(files)
    temp_dir = _make_temp_dir()
    file_list = _make_temp_files(groups, temp_dir)
    main_job = _make_job(groups, file_list, neben_path, primer, amp_size, target)

    if execute:
        job_num = _run_job(main_job)
        results_job = _make_results_job(job_num, target, groups, num_files, species_names, keep)
        _run_job(results_job)
        

def _interpret_primer(primer_id, forward, reverse):
    if primer_id is None:
        return (forward, reverse)

    combos = _unpickle_combos()
    combo = _get_combo_from_id(combos, primer_id)
    return (combo.forward, combo.reverse)


def _unpickle_combos():
    with open("combos.pickle", "rb") as infile:
        combos = pickle.load(infile)
    return combos


def _get_combo_from_id(combos, combo_id):
    for combo in combos:
        if combo.id == combo_id:
            return combo


def _get_neben_path():
    tools_dir = os.path.dirname(os.path.abspath(__file__))
    main_dir = os.path.join(tools_dir, os.pardir)
    pdp_dir = os.path.join(main_dir, "Primer_Design_Pipeline")
    return os.path.join(pdp_dir, "nv2_linux_64")


def _get_files(directory):
    out = []
    for file in glob.glob("{}/**/*.fasta".format(directory), recursive=True):
        out.append(os.path.abspath(file))
    return out


def _split_files(files):
    out = []
    species = {}
    
    time_limit = 600
    group = []
    time_total = 0

    for file in files:
        species_name = os.path.dirname(file).split("/")[-1]
        if species_name in species:
            species[species_name] += 1
        else:
            species[species_name] = 1
        if time_total >= time_limit:
            out.append(group)
            group = []
            time_total = 8

        group.append(file)
        time_total += .5

    if len(group) > 0:
        out.append(group)

    return out, species



def _make_job(groups, file_list, neben_path, primer, amp_size, target):

    outfile_name = "dba_job.sh"

    num_jobs = len(groups)

    with open(outfile_name, "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name=database_amplicon\n")
        outfile.write("#SBATCH --array=0-{}\n".format(num_jobs-1))
        outfile.write("#SBATCH --time=1:00:00\n")
        outfile.write("#SBATCH --mem=10000\n")
        outfile.write("#SBATCH --output=dba_%A.out\n")
        outfile.write("#SBATCH --open-mode=append\n")
        outfile.write("\n")

        outfile.write("FILE_ARRAY=(")
        for file in file_list:
            outfile.write('"{}" '.format(file))
        outfile.write(")\n")
        outfile.write("\n")

        outfile.write("FILE=${FILE_ARRAY[$SLURM_ARRAY_TASK_ID]}\n")
        outfile.write("\n")

        if target:
            species_name = os.path.dirname(groups[0][0]).split("/")[-1]
            subprocess.run('touch temp/{}_output-0.txt\n'.format(species_name), shell=True)
            subprocess.run('touch temp/{}_output-1.txt\n'.format(species_name), shell=True)
            subprocess.run('touch temp/{}_output-2.txt\n'.format(species_name), shell=True)
            
            outfile.write("while read F; do\n")
            outfile.write("\tWC=`{} -m {} -f {} -r {} -g $F | wc -l`".format(
                neben_path,
                amp_size+50,
                primer[0],
                primer[1]))
            outfile.write("\tNUM_AMPS=(${WC// / })\n")
            outfile.write("\tBASENAME=$(basename $F)\n")
            outfile.write("\tif [ $NUM_AMPS -eq 0 ]\n\tthen\n")
            outfile.write("\t\t echo $BASENAME >> temp/{}_output-0.txt\n".format(species_name))
            outfile.write("\telif [ $NUM_AMPS -eq 1 ]\n\tthen\n")
            outfile.write("\t\t echo $BASENAME >> temp/{}_output-1.txt\n".format(species_name))
            outfile.write("\telse\n")
            outfile.write("\t\t echo $BASENAME >> temp/{}_output-2.txt\n".format(species_name))
            outfile.write("\tfi\n")
            

        if not target:
            
            outfile.write("while read F; do\n")
            outfile.write("\tPARENT_ABS=`dirname $F`\n")
            outfile.write("\tPARENT=`basename $PARENT_ABS`\n")
            outfile.write('\tOUTFILE=temp/"$PARENT"_output.txt\n')
            outfile.write("\t{} -m {} -f {} -r {} -g $F | head -1 >> $OUTFILE\n".format(
                neben_path,
                amp_size+50,
                primer[0],
                primer[1]))

        outfile.write("done <$FILE\n")

    return outfile_name


def _make_results_job(job_num, target, groups, num_files, species_names, keep):
    outfile_name = "dbar_job.sh"
    with open(outfile_name, "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name=database_amplicon_results\n")
        outfile.write("#SBATCH --time=20:00\n")
        outfile.write("#SBATCH --mem=2000\n")
        outfile.write("#SBATCH --dependency=afterok:{}\n".format(job_num))
        outfile.write("\n")

        results_name = "results_dba_{}.txt".format(job_num)

        if target:

            num_files = 0
            for group in groups:
                num_files += len(group)

            species_name = os.path.dirname(groups[0][0]).split("/")[-1]

            outfile.write("touch {}\n".format(results_name))

            outfile.write("for i in {0..2}; do\n")
            outfile.write("\tWC=`wc -l temp/{}_output-$i.txt`\n".format(species_name))
            outfile.write("\tNUM_LINES=(${WC// / })\n")
            outfile.write('\tPERCENT=`echo "scale = 4; ($NUM_LINES / {}) * 100" | bc`\n'.format(num_files))
            outfile.write("\tif [ $i -eq 2 ]\n\tthen\n")
            outfile.write('\t\tprintf "$i or more Amplicons\\t$NUM_LINES/{}\\t%.2f%%\\n" $PERCENT >> {}\n'.format(num_files, results_name))
            outfile.write("\telse\n")
            outfile.write('\t\tprintf "$i Amplicons\\t\\t$NUM_LINES/{}\\t%.2f%%\\n" $PERCENT >> {}\n'.format(num_files, results_name))
            outfile.write("\tfi\n")
            outfile.write("done\n")
            outfile.write("\n")

            outfile.write("echo >> {}".format(results_name))
            outfile.write("\n")
            
            outfile.write("for i in {0..2}; do\n")
            outfile.write("\tif [ $i -eq 2 ]\n\tthen\n")
            outfile.write('\t\techo "$i or more Amplicons:" >> {}\n'.format(results_name))
            outfile.write("\telse\n")
            outfile.write('\t\techo "$i Amplicons:" >> {}\n'.format(results_name))
            outfile.write("\tfi\n")
            outfile.write("\tcat temp/{}_output-$i.txt >> {}\n".format(species_name, results_name))
            outfile.write('\techo >> {}\n'.format(results_name))
            outfile.write("done\n")

        if not target:
            output_files = []
            for species in species_names:
                output_files.append((species, "temp/{}_output.txt".format(species)))

            # Make array of species names
            outfile.write("SPECIES=(")
            for file in output_files:
                #outfile.write(file[0] + " ")
                outfile.write('"{}" '.format(file[0]))
            outfile.write(")\n")
            outfile.write("\n")

            # Make array of output file names for each species
            outfile.write("OUTPUT_FILES=(")
            for file in output_files:
                #outfile.write(file[1] + " ")
                outfile.write('"{}" '.format(file[1]))
            outfile.write(")\n")
            outfile.write("\n")

            # Make array of number of files in database for each species
            outfile.write("NUM_FILES=(")
            for file in output_files:
                outfile.write(str(species_names[file[0]]) + " ")
            outfile.write(")\n")
            outfile.write("\n")

            outfile.write("touch {}\n".format(results_name))

            outfile.write("for i in {{0..{}}}; do\n".format(len(output_files)-1))
            outfile.write("\tWC=`wc -l ${OUTPUT_FILES[$i]}`\n")
            outfile.write("\tNUM_LINES=(${WC// / })\n")
            outfile.write("\tif [ $NUM_LINES -gt 0 ]\n\tthen\n")
            outfile.write('\t\tPERCENT=`echo "scale = 4; ($NUM_LINES / ${NUM_FILES[$i]}) * 100" | bc`\n')
            outfile.write('\t\tprintf "${{SPECIES[$i]}}\\t$NUM_LINES/${{NUM_FILES[$i]}}\\t%.2f%%\\n" $PERCENT >> {}\n'.format(results_name))
            outfile.write("\tfi\n")
            outfile.write("done\n")


        outfile.write("if ! [ -s {} ]\nthen\n".format(results_name))
        outfile.write('\techo "No amplicons found" >> {}\n'.format(results_name))
        outfile.write("fi\n")

        if not keep:
            outfile.write("rm -r temp/")

    return outfile_name



def _make_temp_dir():
    tmp_dir = os.path.join(os.getcwd(), "temp")
    if os.path.exists(tmp_dir):
        print("Error - {} directory already exists. Please delete it and try"
              " again".format(tmp_dir))
        sys.exit(1)
    else:
        subprocess.run("mkdir temp", shell=True)

    return tmp_dir


def _make_temp_files(groups, temp_dir):
    file_list = []

    old_dir = os.getcwd()
    os.chdir(temp_dir)
    for idx, group in enumerate(groups):
        file_name = "group-{}.txt".format(idx)
        with open(file_name, "w") as outfile:
            for file in group:
                outfile.write(file + "\n")
        file_list.append(os.path.abspath(file_name))

    os.chdir(old_dir)
    return file_list


def _run_job(job_name):
    out = subprocess.check_output("sbatch {}".format(job_name), shell=True)
    out = out.decode("UTF-8").strip()
    print(out)
    return out.split()[-1]



if __name__ == "__main__":

    # Python paths are annoying
    if __package__ is None:
        sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from Primer_Design_Pipeline import combo_funcs

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="[REQUIRED] Path to database directory")
    parser.add_argument("-e", "--execute", const=True, nargs="?", help="Submit the Slurm job", type=bool, default=False)
    parser.add_argument("-p", "--primer", help="Unique primer ID")
    parser.add_argument("-f", "--forward", help="Forward primer sequence")
    parser.add_argument("-r", "--reverse", help="Reverse primer sequence")
    parser.add_argument("-a", "--amplicon_size", help="Desired amplicon size (default 500)", type=int, default=500)
    parser.add_argument("-t", "--target", const=True, nargs="?", help="Inputted directory is the same species as primer (default False)", type=bool, default=False)
    parser.add_argument("-k", "--keep", const=True, nargs="?", help="Keep temporary files (default False)", type=bool, default=False)

    args = parser.parse_args()

    if args.primer is None and (args.forward is None or args.reverse is None):
        print("Error - must provide a primer ID or a forward and reverse sequence")
        sys.exit(1)

    database_amplicons(args.directory, args.execute, args.primer, args.forward, args.reverse, args.amplicon_size, args.target, args.keep)
