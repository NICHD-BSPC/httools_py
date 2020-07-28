#!/bin/bash
#SBATCH --job-name="httools_py"
#SBATCH --partition="norm"
#SBATCH --time=12:00:00

if
    [ $# -eq 0 ]
then
    echo "ERROR: Config file must be supplied"
    exit 1

fi

if
    [ $HOSTNAME == 'biowulf.nih.gov' ]
then
    # Running on Biowulf
    sbatch --cpus-per-task=4 scripts/WRAPPER_SLURM $1

else
    # Running locally
    scriptdir=$(dirname "$0")
    snakemake -s $scriptdir/scripts/Snakefile --config configfn=$1 -R `snakemake -s $scriptdir/scripts/Snakefile --config configfn=$1 --list-params-changes` --cores=1

fi
