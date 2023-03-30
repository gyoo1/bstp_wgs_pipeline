#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --array=1-12
#SBATCH --mem 20G
#SBATCH -c 8
#SBATCH --time 3:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load fastp

# one way to pass sample list = first passed parameter on command line
# SAMPLELIST=$1

# or, just set the file name directly
SAMPLELIST="file_split${SLURM_ARRAY_TASK_ID}.txt"

# loop through sample list
while IFS=$'\t' read -r -a array; do
fastp -g -w 8 -i ../00-raw/${array[0]}_R1.fastq.gz -I ../00-raw/${array[0]}_R2.fastq.gz -o ./${array[1]}_trimmed_R1.fastq.gz -O ./${array[1]}_trimmed_R2.fastq.gz --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -h ./fastp_reps/${array[1]}.html
done < $SAMPLELIST
