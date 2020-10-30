#!/bin/bash
#SBATCH --job-name="httools"
#SBATCH --partition="norm"
#SBATCH --time=12:00:00

if [ $# -eq 0 ]
  then
    echo "ERROR: Config file must be supplied"
    exit 1

fi

# make logdir
if [[ ! -e logs ]]; then mkdir -p logs; fi

# Run snakemake
(
    time snakemake \
    -s Snakefile \
    -R `snakemake -s Snakefile --config fn=$1 --list-params-changes` \
    -p \
    --directory $PWD \
    -k \
    --rerun-incomplete \
    --jobname "s.{rulename}.{jobid}.sh" \
    -j 999 \
    --cluster-config scripts/clusterconfig.yaml \
    --verbose \
    --cluster 'sbatch {cluster.prefix} --cpus-per-task={threads}  --output=logs/{rule}.o.%j --error=logs/{rule}.e.%j' \
    --use-conda \
    --config fn=$1 \
    --latency-wait=300 \
    ) > "Snakefile.log" 2>&1

SNAKE_PID=$!

finish(){
    echo 'Stopping running snakemake job.'
    kill -SIGINT $SNAKE_PID
    exit 0
}
trap finish SIGTERM

wait $SNAKE_PID

