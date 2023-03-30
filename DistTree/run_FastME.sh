#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=run_fastme
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 20G
#SBATCH -c 4
#SBATCH --time 3:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load fastme

#generate NJ tree
fastme -D 1 -i bstp.dist -o bstp.tree -m b -n b