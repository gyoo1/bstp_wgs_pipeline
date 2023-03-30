#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=angsd_pca
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 20G
#SBATCH -c 8
#SBATCH --time 10:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load python/3.7

source ~/angsd/bin/activate

pcangsd --beagle ./merged_beagles.beagle.gz --out PCAngsd --pcadapt --sites_save --threads 8