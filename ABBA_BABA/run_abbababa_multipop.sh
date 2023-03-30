#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=multipop_Dstats
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 200G
#SBATCH -c 1
#SBATCH --time 24:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load angsd
module load r/4.2.1

#calculate D-statistic
angsd -doAbbababa2 1 -bam bam.filelist -sizeFile sizeFile.size -doCounts 1 -out bam.dstat -anc ../refgenome/LeachsGenome.fasta -ref ../refgenome/LeachsGenome.fasta -useLast 0 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -baq 1 -minQ 20 -minMapQ 20 -p 1

#estimate Z score, p-values
Rscript estAvgError.R angsdFile="bam.dstat" out="Dstat_results_28ind" sizeFile=sizeFile.size errFile=errorList.error nameFile=popNames.name
