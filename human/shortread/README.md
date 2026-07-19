The 1k genomes vcf was downloaded from 

Analysis of the chrX SNP positions present in the [1k genomes project vcf](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.recalibrated_variants.vcf.gz) was performed using closest_het_snp_distance.py, which quantifies the positition of the nearest heterozygous exonic SNP in each individual from the transcriptional start site and transcriptional end site. 

This was then then analysed with the workbook snp_positions2.ipynb, which compares the cummulative number of genes with heterozygous exonic SNPs at different read lengths from the TSS and TES in each individual. This script was also used for calculating the SNP reach per individual using the deduplicated single cell nanopore BAM files from our samples.


