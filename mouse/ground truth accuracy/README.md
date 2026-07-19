Single cell nanopore data was processed using EPI2ME wf-singlecell was filtered for chrX with samtools and deduplicated using UMI-tools. 

scDaisychain was run on the resulting chrX BAM as follows:
```bash
scDaisychain run \
  --bam test_data/deduped.chrX.bam \
  --vcf test_data/CASTxBL6_F1.0p1.vcf.gz \
  --gtf test_data/genes.gtf \
  --outdir scdaisychain_run \
  --original-gene-matrix-dir test_data/gene_raw_feature_bc_matrix \
  --original-transcript-matrix-dir test_data/transcript_raw_feature_bc_matrix \
  --min-reads 10 \
  --lower-cutoff 0.01 \
  --partition-mode weighted \
  --tag-mode count \
  --drop-conflicts \
  --drop-multi-tsv-and-gtf
```

The resulting haplotypes_min10_lc0.01.csv file was used for assessing the accuracy in analysis.ipynb. This script counts the number of SNPs phased correctly vs the total number of SNPs at each stage of the scDaisychain algorithm and overall.

flip_for_groundtruth.py was used to create a corrected groundtruth vcf of all of the SNPs phased by the scDaisychain algorithm.

```bash
python flip_for_groundtruth.py haplotypes_min10_lc0.01.vcf.gz haplotypes_min10_lc0.01.flipped.vcf.gz
```
 scDaisychain was run again from the tag bam stage with the corrected vcf in order to obtain ground truth allelically resolved X chromosome expression.
