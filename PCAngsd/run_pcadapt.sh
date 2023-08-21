#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=pcadapt_pairs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --array=1-5
#SBATCH --mem 20G
#SBATCH -c 8
#SBATCH --time 5:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load python/3.7

pcangsd --beagle ./gl_pair${SLURM_ARRAY_TASK_ID}.beagle.gz --out pair${SLURM_ARRAY_TASK_ID} --pcadapt --sites_save --threads 8