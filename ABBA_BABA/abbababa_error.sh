#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=Dstat_errorfiles
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --array=1-17
#SBATCH --mem 200G
#SBATCH -c 8
#SBATCH --time 5:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load angsd

angsd -doAncError 1 -anc ../refgenome/LeachsGenome.fasta -ref perfectSample.fa -out errorFile_${SLURM_ARRAY_TASK_ID} -bam ../FST/pop_paths/bampath_${SLURM_ARRAY_TASK_ID}.txt