#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=run_ngsDist
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 200G
#SBATCH -c 16
#SBATCH --time 24:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load gsl

#run ngsDist
~/ngs_software/ngsTools/ngsDist/ngsDist \
	--geno ../PCAngsd/merged_beagles.beagle.gz \
	--out bstp.dist --probs --labels ind.label \
	--n_ind 28 --n_sites 15147500 --n_threads 4 \
	--n_boot_rep 50 --boot_block_size 100