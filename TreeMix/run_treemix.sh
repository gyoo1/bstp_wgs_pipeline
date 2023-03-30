#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=treemix
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 200G
#SBATCH -c 8
#SBATCH --time 12:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

for i in {0..5}
do
  treemix -i treemix_input.gz -m $i -o treemix.$i -root ../refgenome/LeachesGenome.fasta -bootstrap -k 500 -noss
done
