import subprocess
import argparse
import pickle
import glob
import sys
import os


def database_amplicons(directory, execute, primer_id, forward, reverse):

    primer = _interpret_primer(primer_id, forward, reverse)
    neben_path = _get_neben_path()
    
    files = _get_files(directory)
    groups = _split_files(files)
    temp_dir = _make_temp_dir()
    file_list = _make_temp_files(groups, temp_dir)
    _make_job(groups, file_list, neben_path, primer)
    

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


# Return dict where key=directory, value=list of lists,
# where each sublist is of files that add up to 20 minutes.
def _split_files(files):
    out = {}
    for key in files:
        out[key] = []
    time_limit = 1200
    
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
            
    return out


def _make_job(groups, file_list, neben_path, primer):

    num_jobs = 0
    for key in groups:
        num_jobs += len(groups[key])

    time = 30 * num_jobs
    hours = time//60
    minutes = time - (hours*60)
    time_str = "{}:{}:00".format(hours, minutes)

    with open("database_neben_job.sh", "w") as outfile:
        outfile.write("#SBATCH --job-name=database_neben\n")
        outfile.write("#SBATCH --array=0-{}\n".format(num_jobs))
        outfile.write("#SBATCH --time={}\n".format(time_str))
        outfile.write("#SBATCH --mem=50000\n")
        outfile.write("\n")

        outfile.write("FILE_ARRAY=(")
        for file in file_list:
            outfile.write(file + " ")
        outfile.write(")\n")
        outfile.write("\n")

        outfile.write("OUTPUT=_output.txt\n")
        outfile.write("FILE=${FILE_ARRAY[$SLURM_ARRAY_TASK_ID]}\n")
        outfile.write("SPLIT_NAME=(${FILE//---/ })\n")
        outfile.write("DIR_NAME=${SPLIT_NAME[0]}\n")
        outfile.write("OUTFILE=$DIR_NAME$OUTPUT\n")
        outfile.write("\n")

        outfile.write("while read f; do\n")
        outfile.write("\t{} -max 500 --primers {}:{} $f >> $OUTFILE\n".format(
            neben_path,
            primer[0],
            primer[1]))
        outfile.write("done <$FILE\n")
        

def _make_temp_dir():
    tmp_dir = os.path.join(os.getcwd(), "tmp")
    if os.path.exists(tmp_dir):
        print("Error - {} directory already exists. Please delete it and try"
              "again".format(tmp_dir))
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
            file_name = "{}---{}".format(directory, idx)
            with open(file_name, "w") as outfile:
                for file in group:
                    outfile.write(file + "\n")
            file_list.append(os.path.abspath(file_name))
            
    os.chdir(old_dir)
    return file_list
            
   
    

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

    args = parser.parse_args()

    if args.primer is None and (args.forward is None or args.reverse is None):
        print("Error - must provide a primer ID or a forward and reverse sequence")
        sys.exit(1)
    
    database_amplicons(args.directory, args.execute, args.primer, args.forward, args.reverse)
