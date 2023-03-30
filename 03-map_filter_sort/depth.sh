#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --array=1-8
#SBATCH --job-name=depth_print
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 20G
#SBATCH -c 8
#SBATCH --time 3:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load samtools

touch depth_${SLURM_ARRAY_TASK_ID}.txt

# Get number of sites, provided they are equal across all files
NSITES=$(awk 'NR==1 {print $NF}' sites.txt)

SAMPLELIST="tr_fl_${SLURM_ARRAY_TASK_ID}.txt"
for SAMPLE in `cat $SAMPLELIST`; do
echo -n "$SAMPLE:" >> depth_${SLURM_ARRAY_TASK_ID}.txt
samtools depth -a ./sorted/${SAMPLE}_sorted_minq20.bam |  awk '{sum+=$3} END { print "Average = ",sum/$NSITES}' >> depth_${SLURM_ARRAY_TASK_ID}.txt
done