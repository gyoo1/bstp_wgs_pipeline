#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --array=1-14
#SBATCH --job-name=angsd_sfs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 200G
#SBATCH -c 1
#SBATCH --time 24:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load angsd

# calculate number of individuals
NIND=`wc -l < pops/pop_${SLURM_ARRAY_TASK_ID}.txt`
MININD=`expr $NIND / 2`

# estimate SAF for each population
angsd -b pop_paths/bampath_${SLURM_ARRAY_TASK_ID}.txt -out 01_angsd_sfs/saf_${SLURM_ARRAY_TASK_ID} \
  -ref ../refgenome/LeachsGenome.fasta -anc ../refgenome/LeachsGenome.fasta \
  -GL 2 -doSaf 1 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -baq 1 \
  -minMapQ 20 -minQ 20 -minInd $MININD