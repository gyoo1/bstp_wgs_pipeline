#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=angsd_2dsfs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 200G
#SBATCH -c 16
#SBATCH --time 24:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load angsd

# calculate 2D SFS for all pairs
for a in {1..17}; do 
	shift
	for b in {1..17}; do 
		realSFS -P 16 01_angsd_sfs/saf_$a.saf.idx 01_angsd_sfs/saf_$b.saf.idx > 02_angsd_sfs/SFS_$a.$b.sfs
	done
done