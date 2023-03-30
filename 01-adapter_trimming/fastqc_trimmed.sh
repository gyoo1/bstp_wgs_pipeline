#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --account=def-vlf
#SBATCH --mem 20G
#SBATCH -c 1
#SBATCH --time 24:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load fastqc

unalias ls
ls ../01-adapter_trimming/*.fastq.gz -1 | xargs -n 1 fastqc -o trimmed_reports
