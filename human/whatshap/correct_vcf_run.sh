#!/bin/bash

pid=$(grep ^Pid /proc/self/status)
corelist=$(grep Cpus_allowed_list: /proc/self/status | awk '{print $2}')
host=$(hostname | sed 's/.gadi.nci.org.au//g')
echo subtask $1 running in $pid using cores $corelist on compute node $host

source ~/.bashrc

sample=$1

daisychain_vcf="/home/913/dk4874/scratch/gdata/scDaisychain_paper/human/output/scDaisychain_modes_multi_filtered_weighted/${sample}/haplotypes_min10_lc0.01.vcf.gz"
comparison="/home/913/dk4874/scratch/gdata/scDaisychain_paper/human/processed_data/whatshap/${sample}/comparison/per_snp_phase_comparison.tsv"
out_vcf="/home/913/dk4874/scratch/gdata/scDaisychain_paper/human/processed_data/whatshap/${sample}/comparison/${sample}.daisychain.whatshap_corrected.vcf.gz"

python correct_vcf.py \
  --daisychain-vcf "$daisychain_vcf" \
  --comparison-tsv "$comparison" \
  --out-vcf "$out_vcf"

bcftools index -t "$out_vcf"

