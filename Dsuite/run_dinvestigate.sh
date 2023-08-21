#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=dinvestigate
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 100G
#SBATCH -c 10
#SBATCH --time 24:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load bcftools

~/Dsuite/Build/Dsuite Dinvestigate -n fd -g -w 2000,400 petrel.vcf.gz sets.txt test_trios.txt
