#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=angsd_fst_win
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 20G
#SBATCH -c 4
#SBATCH --time 12:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load angsd

# calculate FST in 50kb windows with 10kb steps
for c in {1..14}; do
  shift
  for d in {1..14}; do
    realSFS fst stats2 03_angsd_fst/FST_$a.$b.fst.idx -win 50000 -step 10000 -whichFST 1 > 03_angsd_fst/FST_$a.$b.fst.txt -P 4
  done
done