#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --array=1-14
#SBATCH --job-name=tajima_nd
#SBATCH --mem 100G
#SBATCH -c 4
#SBATCH --time 5:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load angsd

# estimate diversity statistics for each scaffold
thetaStat do_stat POP${SLURM_ARRAY_TASK_ID}.thetas.idx

# compute statistics in sliding windows of 10kb with 5kb steps
thetaStat do_stat POP${SLURM_ARRAY_TASK_ID}.thetas.idx -win 10000 -step 5000 -outnames POP${SLURM_ARRAY_TASK_ID}.tajima
