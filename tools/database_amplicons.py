import subprocess
import argparse
import pickle
import glob
import sys
import os


def database_amplicons(directory, execute, primer_id, forward, reverse, amp_size, target):

    primer = _interpret_primer(primer_id, forward, reverse)
    neben_path = _get_neben_path()

    files = _get_files(directory)
    num_files = _get_num_files(files)
    groups = _split_files(files)
    temp_dir = _make_temp_dir()
    file_list = _make_temp_files(groups, temp_dir)
    main_job = _make_job(groups, file_list, neben_path, primer, amp_size, target)

    if execute:
        job_num = _run_job(main_job)
        results_job = _make_results_job(job_num, target, groups, num_files)
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
    return os.path.join(pdp_dir, "neben_linux_64")


def _get_files(directory):
    out = {}
    for file in glob.glob("{}/**/*.fasta".format(directory), recursive=True):
        dir_name = os.path.dirname(file)
        if dir_name not in out:
            out[dir_name] = [os.path.abspath(file)]
        else:
            out[dir_name].append(os.path.abspath(file))
    return out


def _get_num_files(files):
    out = {}
    for directory in files:
        key = directory.split("/")[-1]
        out[key] = len(files[directory])
    return out


# Return dict where key=directory, value=list of lists,
# where each sublist is of files that add up to 20 minutes.
def _split_files(files):
    out = {}
    for key in files:
        out[key] = []
    time_limit = 600

    for key in files:
        group = []
        time_total = 0
        for file in files[key]:
            if time_total >= time_limit:
                out[key].append(group)
                group = []
                time_total = 0

            if "Complete" in file:
                time_total += 3
            elif "Contig" in file:
                time_total += 1
            elif "Scaffold" in file:
                time_total += 1
            else:
                time_total += 3

            group.append(file)
        out[key].append(group)

    return out


def _make_job(groups, file_list, neben_path, primer, amp_size, target):

    outfile_name = "database_amplicon_job.sh"

    num_jobs = 0
    for key in groups:
        num_jobs += len(groups[key])

    '''
    time = 10 * num_jobs
    hours = time//60
    minutes = time - (hours*60)
    time_str = "{}:{}:00".format(hours, minutes)
    '''

    with open(outfile_name, "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name=database_amplicon\n")
        outfile.write("#SBATCH --array=0-{}\n".format(num_jobs-1))
        outfile.write("#SBATCH --time=45:00\n")
        outfile.write("#SBATCH --mem=5000\n")
        outfile.write("\n")

        outfile.write("FILE_ARRAY=(")
        for file in file_list:
            outfile.write('"{}" '.format(file))
        outfile.write(")\n")
        outfile.write("\n")

        outfile.write("OUTPUT=_output\n")
        outfile.write("FILE=${FILE_ARRAY[$SLURM_ARRAY_TASK_ID]}\n")
        outfile.write("SPLIT_NAME=(${FILE//---/ })\n")
        outfile.write("DIR_NAME=${SPLIT_NAME[0]}\n")
        outfile.write("OUTFILE=$DIR_NAME$OUTPUT\n")
        outfile.write("\n")

        if target:
            species_name = list(groups.keys())[0].split("/")[-1]
            subprocess.run('touch tmp/{}_output-0.txt\n'.format(species_name), shell=True)
            subprocess.run('touch tmp/{}_output-1.txt\n'.format(species_name), shell=True)
            subprocess.run('touch tmp/{}_output-2.txt\n'.format(species_name), shell=True)
            
            outfile.write("while read f; do\n")
            outfile.write("\tWC=`{} -max {} --primers {}:{} $f | wc -l`".format(
                neben_path,
                amp_size+50,
                primer[0],
                primer[1]))
            outfile.write("\tNUM_AMPS=(${WC// / })\n")
            outfile.write("\tBASENAME=$(basename $f)\n")
            outfile.write("\tif [ $NUM_AMPS -eq 0 ]\n\tthen\n")
            outfile.write("\t\t echo $BASENAME >> $OUTFILE-0.txt\n")
            outfile.write("\telif [ $NUM_AMPS -eq 1 ]\n\tthen\n")
            outfile.write("\t\t echo $BASENAME >> $OUTFILE-1.txt\n")
            outfile.write("\telse\n")
            outfile.write("\t\t echo $BASENAME >> $OUTFILE-2.txt\n")
            outfile.write("\tfi\n")
            

        if not target:
            
            outfile.write("while read f; do\n")
            outfile.write("\t{} -max {} --primers {}:{} $f | head -1 >> $OUTFILE.txt\n".format(
                neben_path,
                amp_size+50,
                primer[0],
                primer[1]))

        outfile.write("done <$FILE\n")

    return outfile_name


def _make_results_job(job_num, target, groups, num_files):
    outfile_name = "database_amplicon_results_job.sh"
    with open(outfile_name, "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name=database_amplicon_results\n")
        outfile.write("#SBATCH --time=30:00\n")
        outfile.write("#SBATCH --mem=10000\n")
        outfile.write("#SBATCH --dependency=afterok:{}\n".format(job_num))
        outfile.write("\n")

        results_name = "database_amplicon_results.txt"

        if target:

            path = list(groups.keys())[0]
            num_files = 0
            for group in groups[path]:
                num_files += len(group)

            species_name = path.split("/")[-1]
            outfile.write("for i in {0..2}; do\n")
            outfile.write("\tWC=`wc -l tmp/{}_output-$i`".format(species_name))
            outfile.write("\tNUM_LINES=${WC// / })\n")
            outfile.write('\tPERCENT=`echo "scale = 2; ($NUM_LINES / {}) * 100" | bc`\n'.format(num_files))
            outfile.write("\tif [ $i -eq 2 ]\n\tthen\n")
            outfile.write('\t\tprintf "$i or more Amplicons\\t$NUM_LINES/{}\\t$PERCENT%%\\n" >> {}\n'.format(num_files, results_name))
            outfile.write("\telse\n")
            outfile.write('\tprintf "$i Amplicons\\t\\t$NUM_LINES/{}\\t$PERCENT%%\\n" >> {}\n'.format(num_files, results_name))
            outfile.write("\tfi\n")
            outfile.write("done\n")
            outfile.write("\n")
            
            outfile.write("for i in {0..2}; do\n")
            outfile.write("\tif [ $i -eq 2 ]\n\tthen\n")
            outfile.write('\t\techo "$i or more Amplicons:" >> {}\n'.format(results_name))
            outfile.write("\telse\n")
            outfile.write('\t\techo "$i Amplicons:" >> {}\n'.format(results_name))
            outfile.write("\tfi\n")
            outfile.write("\tcat tmp/{}_output-$i.txt >> {}\n".format(species_name, results_name))
            outfile.write('\techo >> {}\n'.format(results_name))
            outfile.write("done\n")

        if not target:
            output_files = []
            for directory in groups:
                species_name = directory.split("/")[-1]
                output_files.append((species_name, "tmp/{}_output.txt".format(species_name)))

            # Make array of species names
            outfile.write("SPECIES=(")
            for file in output_files:
                outfile.write(file[0] + " ")
            outfile.write(")\n")
            outfile.write("\n")

            # Make array of output file names for each species
            outfile.write("OUTPUT_FILES=(")
            for file in output_files:
                outfile.write(file[1] + " ")
            outfile.write(")\n")
            outfile.write("\n")

            # Make array of number of files in database for each species
            outfile.write("NUM_FILES=(")
            for file in output_files:
                outfile.write(str(num_files[file[0]]) + " ")
            outfile.write(")\n")
            outfile.write("\n")

            outfile.write("for i in {{0..{}}}; do\n".format(len(output_files)-1))
            outfile.write("\tWC=`wc -l ${OUTPUT_FILES[$i]}`\n")
            outfile.write("\tNUM_LINES=(${WC// / })\n")
            outfile.write("\tif [ $NUM_LINES -gt 0 ]\n\tthen\n")
            outfile.write('\t\tPERCENT=`echo "scale = 2; ($NUM_LINES / ${NUM_FILES[$i]}) * 100" | bc`\n')
            outfile.write('\t\tprintf "${{SPECIES[$i]}}\\t$NUM_LINES/${{NUM_FILES[$i]}}\\t$PERCENT%%\\n" >> {}\n'.format(results_name))
            outfile.write("\tfi\n")
            outfile.write("done\n")

    return outfile_name



def _make_temp_dir():
    tmp_dir = os.path.join(os.getcwd(), "tmp")
    if os.path.exists(tmp_dir):
        print("Error - {} directory already exists. Please delete it and try"
              " again".format(tmp_dir))
        sys.exit(1)
    else:
        subprocess.run("mkdir tmp", shell=True)

    return tmp_dir


def _make_temp_files(groups, temp_dir):
    file_list = []

    old_dir = os.getcwd()
    os.chdir(temp_dir)
    for directory in groups:
        for idx, group in enumerate(groups[directory]):
            # Weird delimitter "---" that will (hopefully)
            # not show up in any file names.
            file_name = "{}---{}.txt".format(directory.split("/")[-1], idx)
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

    args = parser.parse_args()

    if args.primer is None and (args.forward is None or args.reverse is None):
        print("Error - must provide a primer ID or a forward and reverse sequence")
        sys.exit(1)

    database_amplicons(args.directory, args.execute, args.primer, args.forward, args.reverse, args.amplicon_size, args.target)
