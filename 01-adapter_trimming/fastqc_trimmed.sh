#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --account=def-vlf
#SBATCH --mem 20G
#SBATCH -c 8
#SBATCH --time 24:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load fastqc

unalias ls
ls ./*.fastq.gz -1 | xargs -n 8 fastqc -o ./trimmed_reports