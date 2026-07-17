#!/bin/bash

#PBS -P oo78
#PBS -q normalsr
#PBS -l mem=1000GB
#PBS -l ncpus=208
#PBS -l wd
#PBS -l walltime=48:00:00
#PBS -j oe
#PBS -m ae
#PBS -l storage=gdata/oo78+scratch/oo78+gdata/vn68+gdata/if89
#PBS -l jobfs=800GB

source ~/.bashrc

module load bcftools

awk -F',' '$2 == "female" {print $1}' 1k_sexes.csv > female_samples.txt
# Make a Unix-normalized copy
sed 's/\r$//' 1k_sexes.csv > 1k_sexes.unix.csv

# Recreate the list (case-insensitive, skip header)
awk -F',' 'NR>1 && tolower($2)=="female"{print $1}' 1k_sexes.unix.csv > female_samples.txt

# Sanity check
wc -l female_samples.txt
head female_samples.txt

comm -12 <(sort -u female_samples.txt) \
         <(bcftools query -l /home/913/dk4874/scratch/primateXInactivation/scripts/paper/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.recalibrated_variants.vcf.gz | sort -u) \
> females_in_vcf.txt

wc -l females_in_vcf.txt
head females_in_vcf.txt

bcftools view \
  -S females_in_vcf.txt \
  -Oz \
  -o ALL.chrX.female_only.vcf.gz \
  /home/913/dk4874/scratch/primateXInactivation/scripts/paper/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.recalibrated_variants.vcf.gz

bcftools index ALL.chrX.female_only.vcf.gz
bcftools query -l ALL.chrX.female_only.vcf.gz | wc -l


python closest_het_snp_distances.py \
  --vcf /home/913/dk4874/scratch/primateXInactivation/scripts/paper/ALL.chrX.female_only.vcf.gz \
  --gtf /home/913/dk4874/scratch/primateXInactivation/raw_data/genome/homo_sapiens/genes/genes.gtf \
  --out closest_het_snp_distances.csv \
  --processes 208

