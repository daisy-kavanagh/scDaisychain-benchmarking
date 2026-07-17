#!/bin/bash

pid=$(grep ^Pid /proc/self/status)
corelist=$(grep Cpus_allowed_list: /proc/self/status | awk '{print $2}')
host=$(hostname | sed 's/.gadi.nci.org.au//g')
echo subtask $1 running in $pid using cores $corelist on compute node $host

source ~/.bashrc
sample=$1
ref=="genome/homo_sapiens/fasta/genome.fa"
vcf=/home/913/dk4874/scratch/XCI_RA_aging/processed_data/clair3/${sample}/merge_output_het_biallelic_chrX.vcf.gz
genome_bam=/home/913/dk4874/scratch/vn68_gdata/XCI_RA_genome/processed_data/genome/${sample}/${sample}.reads.pass.aln.coord.bam
outdir=/home/913/dk4874/scratch/gdata/scDaisychain_paper/human/processed_data/whatshap/${sample}/
mkdir -p $outdir


qual20gq10_vcf=/home/913/dk4874/scratch/XCI_RA_aging/processed_data/clair3/${sample}/merge_output_het_biallelic_chrX.QUAL20.gq10.vcf.gz
bcftools view \
  -m2 -M2 \
  -v snps \
  -f PASS \
  -i 'GT="het" && GT~"\|" && FORMAT/DP>=10 && FORMAT/GQ>=10 && QUAL>=20' \
  "$vcf" \
  -Oz -o "$qual20gq10_vcf"

bcftools index -t "$qual20gq10_vcf"


qual20gq10_phased_vcf=${outdir}/merge_output_het_biallelic_chrX.phased.distrust.QUAL20.gq10.vcf.gz


mamba activate whatshap
whatshap phase \
  --reference $ref \
  --chromosome chrX \
  -o $qual20gq10_phased_vcf \
  --only-snvs \
  --ignore-read-groups \
  --distrust-genotypes \
  $qual20gq10_vcf \
  $genome_bam \
  
bcftools index -t $qual20gq10_phased_vcf

daisychain_vcf=/home/913/dk4874/scratch/gdata/scDaisychain_paper/human/output/scDaisychain_modes_multi_filtered_weighted/${sample}/haplotypes_min10_lc0.01.vcf.gz
comparison_out=/home/913/dk4874/scratch/gdata/scDaisychain_paper/human/processed_data/whatshap/${sample}/comparison_qual20gq10
mkdir -p $comparison_out
python compare_phase_concordance.py --vcf1 $qual20gq10_phased_vcf --vcf2 $daisychain_vcf --outdir $comparison_out



