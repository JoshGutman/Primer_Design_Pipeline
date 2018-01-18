#!/bin/bash
#SBATCH --job-name=primer_design_pipeline
#SBATCH -c 4
#SBATCH --time=3:00:00
#SBATCH --mem=50000

python ../../primer_design_pipeline.py -t multifasta.fasta -d ../ -c primer3_config.txt 
-g target_list.txt -u 325 -l 275 -i 98 -ol .2 -na 0 -mg 3
