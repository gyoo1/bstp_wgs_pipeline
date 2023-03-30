#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=allele_freq_bam
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --array=1-6
#SBATCH --mem 200G
#SBATCH -c 1
#SBATCH --time 12:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load angsd

# calculate number of individuals
NIND=`wc -l < pops/pop_${SLURM_ARRAY_TASK_ID}.txt`
MININD=`expr $NIND / 2`

angsd -bam pops/pop_${SLURM_ARRAY_TASK_ID}.txt -GL 2 -doMaf 3 -doMajorMinor 1 \
	-doCounts 1 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd $MININD \
	-out freq_pop_${SLURM_ARRAY_TASK_ID} -sites sites.txt -fai ../refgenome/LeachsGenome.fasta.fai
