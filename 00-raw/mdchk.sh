#!/bin/bash
#SBATCH --job-name=hashcheck
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --account=def-vlf
#SBATCH --mem 4G
#SBATCH -c 1
#SBATCH --time 1:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

# remove any unusual ls behaviours

unalias ls

echo "Number of files: " > check.txt
ls -la | wc -l >> check.txt

# loop through all checksums
for hash in $(ls *.md5 -1); do
  md5sum -c $hash >> check.txt
done