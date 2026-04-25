#!/bin/bash
#
#SBATCH --job-name=ParaChrSNP
#SBATCH -c 15 # Number of cores
#SBATCH -N 1 # Run the Snakemake controller on one node
#SBATCH --mem=20G # Memory pool for all cores in MB (see also --mem-per-cpu)
#SBATCH -p pNormal # Partition to submit to
#SBATCH --qos=normal
#SBATCH -o snake.out # File to which STDOUT will be written
#SBATCH -e snake.err # File to which STDERR will be written

cd "$(dirname "$0")"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate parachrsnp

snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --jobs 5 \
    --keep-going \
    --cluster-config cluster.yaml \
    --cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n} -c {cluster.c} --mem-per-cpu {cluster.mempercpu}"
