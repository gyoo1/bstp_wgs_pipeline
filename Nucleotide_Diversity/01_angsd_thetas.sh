#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --array=1-14
#SBATCH --job-name=angsd_thetas
#SBATCH --mem 200G
#SBATCH -c 8
#SBATCH --time 5:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load angsd

# compute thetas for each site
realSFS saf2theta ../FST/01_angsd_sfs/saf_${SLURM_ARRAY_TASK_ID}.saf.idx -sfs ../FST/01_angsd_sfs/SFS_${SLURM_ARRAY_TASK_ID} -outname POP${SLURM_ARRAY_TASK_ID} -fold 1 -P 8
