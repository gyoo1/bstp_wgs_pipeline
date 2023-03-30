#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --array=1-8
#SBATCH --job-name=align_rt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 20G
#SBATCH -c 8
#SBATCH --time 8:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load bwa
module load samblaster
module load samtools

# $SLURM_ARRAY_TASK_ID is the id of the current script copy

SAMPLELIST="tr_fl_${SLURM_ARRAY_TASK_ID}.txt"

for SAMPLE in `cat $SAMPLELIST`; do

echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

bwa mem -t 8 ../refgenome/LeachesGenome.fasta ../01-adapter_trimming/${SAMPLE}_R1.fastq.gz ../01-adapter_trimming/${SAMPLE}_R2.fastq.gz | samblaster --removeDups | samtools view -h -b -@8 -o ./${SAMPLE}_aligned.bam

done
