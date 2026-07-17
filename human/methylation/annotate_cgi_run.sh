#!/bin/bash
set -euo pipefail

sample="${1:?Usage: $0 SAMPLE_NAME}"

module load bedtools/2.31.0

base="/home/913/dk4874/scratch/vn68_gdata/XCI_RA_genome/processed_data/bedmethyl"
indir="${base}/${sample}"

in="${indir}/${sample}.reads.pass.modkit.coordA.cpg.chrX.bed"
in_gz="${indir}/${sample}.reads.pass.modkit.coordA.cpg.sorted.bed.gz"
out="${indir}/${sample}.reads.pass.modkit.coordA.cpg.chrX_annotated_with_CGI_shores_TSS.bed"

cpg_islands="/home/913/dk4874/scratch/gdata/genomes/cellranger/homo_sapiens/refdata-gex-GRCh38-2024-A/genes/chrX_cpg_islands.bed"
tss="tss="/home/913/dk4874/scratch/X_Inactivation/raw_data/genome/homo_sapiens_cellranger/genes/promoters_TSS_1kb_gene.bed""

[[ -d "$indir" ]] || { echo "ERROR: missing sample dir: $indir" >&2; exit 1; }
[[ -s "$cpg_islands" ]] || { echo "ERROR: missing CpG islands: $cpg_islands" >&2; exit 1; }
[[ -s "$tss" ]] || { echo "ERROR: missing TSS BED: $tss" >&2; exit 1; }

if [[ ! -s "$in" ]]; then
  [[ -s "$in_gz" ]] || { echo "ERROR: missing input: $in or $in_gz" >&2; exit 1; }
  echo "Creating chrX-only BED: $in"
  zcat "$in_gz" | awk 'BEGIN{OFS="\t"} $1=="chrX"' > "$in"
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

echo "Input:  $in"
echo "Output: $out"

# CpG shores = +/- 2 kb around CGI, excluding CGI itself
awk 'BEGIN{OFS="\t"} {
  s=$2-2000; if (s<0) s=0
  print $1,s,$3+2000,"CGI_shore"
}' "$cpg_islands" \
  | bedtools subtract -a - -b "$cpg_islands" \
  > "$tmp/cpg_shores.bed"

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"CGI"}' "$cpg_islands" \
  > "$tmp/cpg_islands.named.bed"

cat "$tmp/cpg_islands.named.bed" "$tmp/cpg_shores.bed" \
  | sort -k1,1 -k2,2n \
  > "$tmp/cgi_and_shores.bed"

# Add row ID as col 19 to preserve duplicate h/m rows per CpG
awk 'BEGIN{OFS="\t"} {print $0, NR}' "$in" > "$tmp/input.with_id.bed"

# Add CGI/shore annotation as col 19, keeping original 18 columns
bedtools intersect -a "$tmp/input.with_id.bed" -b "$tmp/cgi_and_shores.bed" -loj \
  | awk 'BEGIN{OFS="\t"}
    {
      id=$19
      annot = ($(NF-3)=="." ? "non-CGI" : $NF)

      if (!(id in row)) {
        row[id]=$0
        cgi[id]=annot
        order[++n]=id
      } else if (annot=="CGI") {
        cgi[id]="CGI"
      } else if (cgi[id]=="non-CGI" && annot!="non-CGI") {
        cgi[id]=annot
      }
    }
    END{
      for (i=1;i<=n;i++) {
        id=order[i]
        split(row[id],a,OFS)
        for (j=1;j<=18;j++) printf "%s%s", a[j], OFS
        print cgi[id], id
      }
    }' > "$tmp/with_cgi.bed"

# Add TSS/promoter annotation as final column, then remove row ID
bedtools intersect -a "$tmp/with_cgi.bed" -b "$tss" -loj \
  | awk 'BEGIN{OFS="\t"}
    {
      id=$20

      # A has 20 cols. B is BED6, so B fields are 21-26.
      # B name is col 24.
      tss_name = ($21=="." ? "." : $24)

      if (!(id in row)) {
        row[id]=$0
        tss[id]=tss_name
        order[++n]=id
      } else if (tss[id]=="." && tss_name!=".") {
        tss[id]=tss_name
      } else if (tss_name!="." && index("," tss[id] ",", "," tss_name ",")==0) {
        tss[id]=tss[id] "," tss_name
      }
    }
    END{
      for (i=1;i<=n;i++) {
        id=order[i]
        split(row[id],a,OFS)
        for (j=1;j<=19;j++) printf "%s%s", a[j], OFS
        print tss[id]
      }
    }' > "$out"

echo "Done."
wc -l "$in" "$out"
head "$out"