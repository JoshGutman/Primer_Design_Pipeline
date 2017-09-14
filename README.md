# Primer_Design_Pipeline

## Requirements

* Python 3.5 or newer
* BioPython
* blastn
* muscle

## Usage manual

### Flags
```
python primer_design_pipeline.py [flags]

Required:
    -t str
        Path to target multifasta
    -d str
        Path to reference directory
    -c str
        Path to primer3 config file
    -g str
        Path to genome target list

Optional:
    -r str
        Path to amplicon reference assembly (combined.seqs)
    -l int
        Lower bound of amplicon size (150)
    -u int
        Upper bound of amplicon size (250)
    -i float
        Threshold percentage to consider SNP a degeneracy (98.0)
    -ol float
        Oligo concentration (uM) for calculating Tm (0.25)
    -na float
        Na+ concentration (mM) for calculating Tm (50.0)
    -mg float
        Mg++ concentration (mM) for calculating Tm (0.00)
    -k bool
        Keep all temporary files (False)
```

### Instructions
1. Copy or link all reference fastas into a directory (to be referred to as the reference directory).
2. If you don't have a target list, run `python tools/make_target_list.py -d /path/to/reference/directory`.
3. Copy or link all outgroups to reference directory.
4. Make primer3 config file (a template can be found in `tools/`).
5. Make target multifasta.
6. Move target list, primer3 config file, and multifasta to their own directory.
7. Make sure that the only files in the reference directory are the reference genomes and the outgroups.
8. Load appropriate modules (`anaconda/3.latest`, `blast+/2.2.29`, `muscle/3.8.31`).
9. Run `python primer_design_pipeline.py [flags]` on command-line or as a SLURM job (recommended).

### Expected job stats
Specifications:
```
Fasta files in reference directory: 540
Targets in target multifasta: 2
Processors used: 4
```

Stats:
```
Maximum memory usage: 693 MB
Elapsed Time: 00:06:57
```

### Output
In the directory that you ran the job from, there will be a file called best_primers.txt. This file contains an overview of the three best primers from each target in the target multifasta.

In the directory that you ran the job from, there will be a sub-directory for each target in the target multifasta. In each sub-directory, there will be four files (if the -k flag is False):
1. candidate_primers.txt -- Same as best_primers.txt, except with all of the possible primers instead of the three best ones.
2. ordering_info.txt -- Ordering information for each possible primer.
3. primer_conflicts.txt -- Quality check that shows which primers conflict with each other and at what lengths.
4. primers.blast.out -- Quality check that BLASTs all primers against the combined BLAST database (combined.seqs).

The slurm outfile will include a lot of junk at the beginning -- this is due to neben. If an error occurs, it will be at the bottom of the outfile.