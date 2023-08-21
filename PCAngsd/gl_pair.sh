#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=gl_pairs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --array=1-6
#SBATCH --mem 200G
#SBATCH -c 1
#SBATCH --time 24:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load angsd

# calculate number of individuals
NIND=$(wc -l pairs/pair${SLURM_ARRAY_TASK_ID}.txt)
MININD=$(echo "($NIND / 2)" | bc)
MAXDEPTH=$(echo "($NIND * 8.7 / 1)" | bc)

# generate beagle files for PCA
angsd -b pairs/pair${SLURM_ARRAY_TASK_ID}.txt \
  -out gl_pair${SLURM_ARRAY_TASK_ID} -ref ../refgenome/LeachsGenome.fasta \
  -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doCounts 1 -doCheck 0 -nThreads 1 \
  -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
  -minInd $MININD -setMinDepth $NIND -setMaxDepth $MAXDEPTH -SNP_pval 1e-6 -minMaf 0.01