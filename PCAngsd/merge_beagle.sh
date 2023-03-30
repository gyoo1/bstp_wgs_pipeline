#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=merge_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 4G
#SBATCH -c 1
#SBATCH --time 6:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

# your beagle file list, omitting the first file (which will retain its header)
SAMPLELIST="beagle_files.txt"
HFILE="./angsd_output_reg1/reg1_angsd.beagle"

# generate tmp file with header
zcat ${HFILE}.gz > ${HFILE}.tmp

# open and clear list file for tmps, print first file name to list
printf "${HFILE}.tmp\n" > beagle_tmps.list

# loop through all remaining files to be merged, remove the header, then save to tmp
# print tmp file name to list file
for SAMPLE in `cat $SAMPLELIST`; do
  zcat ${SAMPLE}.gz | tail -n +2 > ${SAMPLE}.tmp
  printf "${SAMPLE}.tmp\n" >> beagle_tmps.list
done

# use xargs to combine all tmp files on the list
xargs cat < beagle_tmps.list > merged_beagles.beagle

# gzip the merged file
gzip merged_beagles.beagle
