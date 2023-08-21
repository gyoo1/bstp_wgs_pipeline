#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=angsd_fst
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 200G
#SBATCH -c 16
#SBATCH --time 24:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load angsd

# calculate FST
for a in {1..14}; do
  shift
  for b in {1..14}; do
    realSFS fst index 01_angsd_sfs/saf_$a.saf.idx 01_angsd_sfs/saf_$b.saf.idx \
      -sfs 02_angsd_sfs/SFS_$a.$b.sfs -fstout 03_angsd_fst/FST_$a.$b -whichFST 1 -P 16 
  done
done