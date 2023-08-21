#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=dtrios_reg1-7
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 100G
#SBATCH -c 10
#SBATCH --time 30:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load bcftools

~/Dsuite/Build/Dsuite Dtrios -c -n dstats -g -t constree.newick reg1-7.vcf.gz sets.txt
