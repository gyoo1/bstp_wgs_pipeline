#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --array=1-14
#SBATCH --job-name=angsd_sfs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 200G
#SBATCH -c 12
#SBATCH --time 12:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load angsd

# estimate SFS
realSFS -P 12 01_angsd_sfs/saf_${SLURM_ARRAY_TASK_ID}.saf.idx > 01_angsd_sfs/SFS_${SLURM_ARRAY_TASK_ID}