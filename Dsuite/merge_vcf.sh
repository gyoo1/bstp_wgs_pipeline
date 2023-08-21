#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=merge_vcfs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=20gy11@queensu.ca
#SBATCH --mem 50G
#SBATCH -c 1
#SBATCH --time 24:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load bcftools

bcftools view angsd_output_reg2/reg2_angsd.bcf | tail -n+1743 > reg2_vcf.txt
bcftools view angsd_output_reg3/reg3_angsd.bcf | tail -n+1743 > reg3_vcf.txt
bcftools view angsd_output_reg4/reg4_angsd.bcf | tail -n+1743 > reg4_vcf.txt
bcftools view angsd_output_reg5/reg5_angsd.bcf | tail -n+1743 > reg5_vcf.txt
bcftools view angsd_output_reg6/reg6_angsd.bcf | tail -n+1743 > reg6_vcf.txt
bcftools view angsd_output_reg7/reg7_angsd.bcf | tail -n+1743 > reg7_vcf.txt
bcftools view angsd_output_reg8/reg8_angsd.bcf | tail -n+1743 > reg8_vcf.txt
bcftools view angsd_output_reg9/reg9_angsd.bcf | tail -n+1743 > reg9_vcf.txt
bcftools view angsd_output_reg10/reg10_angsd.bcf | tail -n+1743 > reg10_vcf.txt

cat vcf_header.txt reg1_vcf.txt reg2_vcf.txt reg3_vcf.txt reg4_vcf.txt reg5_vcf.txt reg6_vcf.txt reg7_vcf.txt reg8_vcf.txt reg9_vcf.txt reg10_vcf.txt | gzip > petrel.vcf.gz
