#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=angsd_grid
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --array=1-10
#SBATCH --mem 200G
#SBATCH -c 1
#SBATCH --time 24:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load angsd

SL="../03-map_filter_sort/sorted/region${SLURM_ARRAY_TASK_ID}.rf"
FOLD="angsd_output_reg${SLURM_ARRAY_TASK_ID}"

# make directory for angsd output
mkdir $FOLD

# generate beagle files for PCA
angsd -b bampath.txt -out "${FOLD}/reg${SLURM_ARRAY_TASK_ID}_angsd" \
	-ref ../refgenome/LeachsGenome.fasta -rf $SL \
	-GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doCounts 1 -doCheck 0 -nThreads 1 \
	-uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 14 -setMinDepth 14 -setMaxDepth 190 \
	-SNP_pval 1e-6 -minMaf 0.01

