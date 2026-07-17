#!/usr/bin/env python3
import argparse, gzip, csv, sys
from collections import defaultdict, namedtuple
from multiprocessing import Pool, cpu_count
from typing import List, Tuple, Dict, Optional
import pysam

Gene = namedtuple("Gene", ["gene_id", "chrom", "strand", "start", "end", "exons_merged", "tss", "tes", "exonic_length"])

# ---------------------- GTF parsing utils ----------------------
def parse_attrs(attr: str) -> Dict[str, str]:
    out = {}
    for kv in attr.strip().split(";"):
        kv = kv.strip()
        if not kv:
            continue
        if " " in kv:
            k, v = kv.split(" ", 1)
            out[k] = v.strip().strip('"')
    return out

def parse_gtf_chrX(gtf_path: str, wanted_chroms=("X","chrX","23")) -> List[Gene]:
    """
    Parse only chrX genes; build merged exons, TSS/TES from gene spans.
    Returns Gene tuples with 0-based half-open coordinates.
    """
    genes_basic: Dict[str, Dict] = {}
    gene_exons: Dict[str, Dict[str, List[Tuple[int,int]]]] = defaultdict(lambda: defaultdict(list))
    tx2gene: Dict[str,str] = {}

    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            if chrom not in wanted_chroms:
                continue
            start_i = int(start) - 1  # 0-based half-open
            end_i   = int(end)
            a = parse_attrs(attrs)

            if feature == "gene":
                gene_id = a.get("gene_id")
                if not gene_id:
                    continue
                genes_basic[gene_id] = {"chrom": chrom, "strand": strand, "start": start_i, "end": end_i}

            elif feature == "transcript":
                gene_id = a.get("gene_id"); tx_id = a.get("transcript_id")
                if gene_id and tx_id:
                    tx2gene[tx_id] = gene_id

            elif feature == "exon":
                tx_id = a.get("transcript_id")
                gene_id = a.get("gene_id") or (tx2gene.get(tx_id) if tx_id else None)
                if not gene_id:
                    continue
                gene_exons[gene_id]["chrom"]  = chrom
                gene_exons[gene_id]["strand"] = gene_exons[gene_id].get("strand", strand)
                gene_exons[gene_id].setdefault("exons", []).append((start_i, end_i))

    genes: List[Gene] = []
    for gene_id, basic in genes_basic.items():
        chrom  = basic["chrom"]
        strand = basic["strand"]
        gstart = basic["start"]
        gend   = basic["end"]

        exons = gene_exons.get(gene_id, {}).get("exons", [])
        exons_merged: List[Tuple[int,int]] = []
        if exons:
            exons_sorted = sorted(exons)
            cs, ce = exons_sorted[0]
            for s,e in exons_sorted[1:]:
                if s <= ce:
                    ce = max(ce, e)
                else:
                    exons_merged.append((cs, ce))
                    cs, ce = s, e
            exons_merged.append((cs, ce))
            # clip defensively to gene span
            exons_merged = [(max(gstart,s), min(gend,e)) for s,e in exons_merged if max(gstart,s) < min(gend,e)]

        if strand == "+":
            tss = gstart
            tes = gend - 1
        else:
            tss = gend - 1
            tes = gstart

        exonic_length = sum(e - s for s, e in exons_merged)
        genes.append(Gene(gene_id, chrom, strand, gstart, gend, exons_merged, tss, tes, exonic_length))
    return genes

# ---------------------- coordinate helpers ----------------------
def build_exonic_index(gene: Gene) -> Tuple[List[Tuple[int,int,int]], int]:
    """
    Returns (segments, total_length) where segments is a list of (start, end, cum_from_tss)
    ordered from TSS to TES in transcription order. All positions are 0-based half-open.
    """
    exons = gene.exons_merged
    if not exons:
        return [], 0
    exons_ord = exons if gene.strand == "+" else list(reversed(exons))
    segments = []
    cum = 0
    for s,e in exons_ord:
        segments.append((s,e,cum))
        cum += (e - s)
    return segments, cum

def exonic_offset_from_tss(gene: Gene, segments, pos: int) -> Optional[int]:
    """
    If pos lies in an exon, return its offset (0-based) from the TSS along exonic sequence.
    Otherwise return None.
    """
    if not segments:
        return None
    if gene.strand == "+":
        for s,e,cum in segments:
            if s <= pos < e:
                return cum + (pos - s)
    else:
        for s,e,cum in segments:
            if s <= pos < e:
                return cum + ((e - 1) - pos)
    return None

# ---------------------- VCF helpers ----------------------
def pick_chrom_name_for_X(vcf: pysam.VariantFile) -> Optional[str]:
    """
    Return the contig name to use for chrX based on the VCF header.
    Common cases: 'X', 'chrX', '23'. Otherwise None.
    """
    contigs = set(vcf.header.contigs.keys())
    if "X" in contigs:
        return "X"
    if "chrX" in contigs:
        return "chrX"
    if "23" in contigs:
        return "23"
    return None

def vcf_has_samples(vcf_path: str) -> bool:
    with pysam.VariantFile(vcf_path) as vf:
        return len(vf.header.samples) > 0

def confirm_fetchable_x(vcf_path: str, x_alias: str) -> None:
    """
    Ensure that fetching on x_alias yields at least something or is fetchable.
    Raises SystemExit with a helpful message if not.
    """
    with pysam.VariantFile(vcf_path) as vf:
        try:
            it = vf.fetch(x_alias)
        except ValueError:
            first20 = list(vf.header.contigs.keys())[:20]
            raise SystemExit(
                f"chr alias '{x_alias}' not fetchable in VCF. "
                f"First 20 contigs in header: {first20}"
            )
        # just attempt to iterate a few records to confirm presence
        probe = 0
        for _ in it:
            probe += 1
            if probe > 5:
                break
        # It's fine if probe==0 (e.g., sparse regions will still fetch later by gene windows).

# ---------------------- genotype logic ----------------------
def is_informative_on_X(gt, treat_haploid_alt_as_informative: bool) -> bool:
    """
    Returns True for:
      - diploid heterozygous calls (len == 2 and alleles differ)
      - (optional) haploid ALT calls on X when treat_haploid_alt_as_informative=True
    """
    if gt is None:
        return False
    a = [x for x in gt if x is not None]
    if len(a) == 2:
        return a[0] != a[1]
    if len(a) == 1 and treat_haploid_alt_as_informative:
        return a[0] == 1
    return False

def is_snp_record(rec) -> bool:
    if rec.ref is None or len(rec.ref) != 1:
        return False
    if not rec.alts:
        return False
    return all(len(a) == 1 for a in rec.alts)

# ---------------------- worker ----------------------
def distances_for_sample(args):
    sample, vcf_path, genes, x_name_map, haploid_alt_informative = args
    rows = []
    informative_genes = 0

    vcf = pysam.VariantFile(vcf_path)
    if sample not in vcf.header.samples:
        vcf.close()
        return rows, sample, informative_genes, len(genes), 0, 0

    # precompute exonic segments per gene
    gene_segments = {}
    for g in genes:
        segs, total_len = build_exonic_index(g)
        gene_segments[g.gene_id] = (segs, total_len)

    # ploidy probe (lightweight)
    dip_gt = hap_gt = 0
    try:
        for rec in vcf.fetch(x_name_map.get("X", "X")):
            sm = rec.samples.get(sample)
            if sm is None:
                continue
            gt = sm.get("GT")
            if gt is None:
                continue
            a = [x for x in gt if x is not None]
            if len(a) == 1:
                hap_gt += 1
            elif len(a) == 2:
                dip_gt += 1
            if dip_gt + hap_gt >= 5000:
                break
    except Exception:
        # ignore probe errors; not essential
        pass

    for g in genes:
        # harmonize chr name to VCF
        vcf_chr = x_name_map.get(g.chrom, g.chrom)
        try:
            it = vcf.fetch(vcf_chr, g.start, g.end)
        except ValueError:
            # chrX alias mismatch for this label; skip
            continue

        best_tss_dna = best_tes_dna = None
        best_tss_rna = best_tes_rna = None
        segs, total_ex_len = gene_segments[g.gene_id]

        saw_informative_here = False

        for rec in it:
            if not is_snp_record(rec):
                continue
            sm = rec.samples.get(sample)
            if sm is None:
                continue
            if not is_informative_on_X(sm.get("GT"), haploid_alt_informative):
                continue

            pos = rec.pos - 1  # 0-based
            # DNA distances (genomic)
            d_tss = abs(pos - g.tss)
            d_tes = abs(pos - g.tes)
            if best_tss_dna is None or d_tss < best_tss_dna:
                best_tss_dna = d_tss
            if best_tes_dna is None or d_tes < best_tes_dna:
                best_tes_dna = d_tes

            # RNA distances (exonic only)
            if segs:
                off = exonic_offset_from_tss(g, segs, pos)
                if off is not None:
                    d_tss_rna = off
                    d_tes_rna = max(0, total_ex_len - 1 - off)
                    if best_tss_rna is None or d_tss_rna < best_tss_rna:
                        best_tss_rna = d_tss_rna
                    if best_tes_rna is None or d_tes_rna < best_tes_rna:
                        best_tes_rna = d_tes_rna

            saw_informative_here = True

        def fmt(x): return "" if x is None else int(x)
        rows.append((sample, g.gene_id, fmt(best_tss_rna), fmt(best_tes_rna), fmt(best_tss_dna), fmt(best_tes_dna)))
        if saw_informative_here:
            informative_genes += 1

    vcf.close()
    return rows, sample, informative_genes, len(genes), dip_gt, hap_gt

# ---------------------- main ----------------------
def main():
    ap = argparse.ArgumentParser(
        description="chrX-only: closest heterozygous (or informative) SNP distances to TSS/TES per sample×gene "
                    "(RNA exon-only and DNA genomic)."
    )
    ap.add_argument("--vcf", required=True, help="BGZipped + indexed VCF with genotypes (chrX present)")
    ap.add_argument("--gtf", required=True, help="GTF with gene/exon features")
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument("--processes", type=int, default=max(1, cpu_count() // 2))
    ap.add_argument("--haploid-alt-informative", action="store_true",
                    help="Count haploid ALT calls on chrX as informative (useful for male samples).")
    args = ap.parse_args()

    # GTF -> chrX genes
    genes = parse_gtf_chrX(args.gtf, wanted_chroms=("X", "chrX", "23"))
    if not genes:
        raise SystemExit("No chrX genes parsed from GTF. Check chromosome labels (X/chrX/23).")
    print(f"[info] Parsed {len(genes)} chrX genes from GTF.", file=sys.stderr)

    # VCF sanity: must have samples (not sites-only)
    if not vcf_has_samples(args.vcf):
        raise SystemExit("VCF has no sample columns (sites-only VCF). Provide a genotyped VCF.")

    # Determine chrX alias and confirm fetchability
    with pysam.VariantFile(args.vcf) as vf:
        x_alias = pick_chrom_name_for_X(vf)
        if x_alias is None:
            first20 = list(vf.header.contigs.keys())[:20]
            raise SystemExit(
                "Could not determine chrX contig name in VCF header "
                f"(tried 'X', 'chrX', '23'). First 20 contigs: {first20}"
            )
    confirm_fetchable_x(args.vcf, x_alias)
    print(f"[info] Using chrX alias '{x_alias}' from VCF header.", file=sys.stderr)

    # Build mapping from GTF chrom labels to VCF alias
    x_name_map = {"X": x_alias, "chrX": x_alias, "23": x_alias}

    # Samples
    with pysam.VariantFile(args.vcf) as vcf2:
        samples = list(vcf2.header.samples)
    if not samples:
        raise SystemExit("No samples in VCF after header parse.")

    print(f"[info] Found {len(samples)} samples in VCF.", file=sys.stderr)

    # Tasks
    tasks = [
        (s, args.vcf, genes, x_name_map, args.haploid_alt_informative)
        for s in samples
    ]

    # Run
    all_rows = []
    summaries = []
    with Pool(processes=args.processes) as pool:
        for rows, sample, informative_genes, n_genes, dip_gt, hap_gt in pool.imap_unordered(distances_for_sample, tasks, chunksize=1):
            all_rows.extend(rows)
            summaries.append((sample, informative_genes, n_genes, dip_gt, hap_gt))
            print(f"[summary] {sample}: informative genes {informative_genes}/{n_genes}; "
                  f"X ploidy probe (diploid GTs ~{dip_gt}, haploid GTs ~{hap_gt})",
                  file=sys.stderr)

    # Write CSV
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "gene", "TSS_distance_RNA", "TES_distance_RNA", "TSS_distance_DNA", "TES_distance_DNA"])
        w.writerows(all_rows)

    print(f"[done] Wrote {len(all_rows)} sample×gene rows to {args.out} across {len(samples)} samples and {len(genes)} chrX genes.",
          file=sys.stderr)

if __name__ == "__main__":
    main()
