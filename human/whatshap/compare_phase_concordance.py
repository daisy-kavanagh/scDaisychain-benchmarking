#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import argparse
import math
import pandas as pd
import pysam
import sys

DEFAULT_VCF1 = sys.argv[1]

DEFAULT_VCF2 = sys.argv[2]


def choose_sample(vcf, sample=None):
    samples = list(vcf.header.samples)
    if sample is not None:
        if sample not in samples:
            raise ValueError(f"Sample {sample!r} not found. Available samples: {samples}")
        return sample
    if len(samples) == 0:
        raise ValueError("VCF has no samples.")
    if len(samples) > 1:
        print(f"[WARN] Multiple samples found; using first sample: {samples[0]}")
    return samples[0]


def is_biallelic_snp(rec):
    return (
        rec.ref is not None
        and rec.alts is not None
        and len(rec.alts) == 1
        and len(rec.ref) == 1
        and len(rec.alts[0]) == 1
    )


def get_phase_code(call):
    """
    Return phase code for phased het biallelic genotype.

    0|1 -> 0
    1|0 -> 1

    Returns None for missing, unphased, homozygous, multiallelic, etc.
    """
    gt = call.get("GT")

    if gt is None or len(gt) != 2:
        return None

    if None in gt:
        return None

    if not getattr(call, "phased", False):
        return None

    if gt == (0, 1):
        return 0
    if gt == (1, 0):
        return 1

    return None


def get_ps(call):
    """
    Extract PS field if present.
    """
    try:
        ps = call.get("PS")
    except Exception:
        return None

    if ps is None:
        return None

    if isinstance(ps, (tuple, list)):
        if len(ps) == 0:
            return None
        ps = ps[0]

    return str(ps)


def load_vcf2_sites(vcf_path, sample=None, snps_only=True):
    """
    Load second VCF as a lookup table by exact variant identity:
    chrom, pos, ref, alt.
    """
    vcf = pysam.VariantFile(str(vcf_path))
    sample = choose_sample(vcf, sample)

    sites = {}

    n_total = 0
    n_kept = 0

    for rec in vcf.fetch():
        n_total += 1

        if snps_only and not is_biallelic_snp(rec):
            continue

        call = rec.samples[sample]
        code = get_phase_code(call)

        if code is None:
            continue

        key = (rec.contig, rec.pos, rec.ref, rec.alts[0])
        sites[key] = {
            "chrom": rec.contig,
            "pos": rec.pos,
            "ref": rec.ref,
            "alt": rec.alts[0],
            "code": code,
            "gt": "|".join(map(str, call.get("GT"))),
        }

        n_kept += 1

    print(f"[VCF2] sample={sample}")
    print(f"[VCF2] records read: {n_total:,}")
    print(f"[VCF2] phased het SNPs kept: {n_kept:,}")

    return sites, sample


def compare_phase_sets(vcf1_path, vcf2_sites, sample=None, snps_only=True):
    """
    Use VCF1 phase sets as reference.
    Compare shared phased het SNPs to VCF2.
    """
    vcf1 = pysam.VariantFile(str(vcf1_path))
    sample = choose_sample(vcf1, sample)

    phase_sets = defaultdict(list)

    n_total = 0
    n_vcf1_phased_het = 0
    n_no_ps = 0
    n_shared = 0

    for rec in vcf1.fetch():
        n_total += 1

        if snps_only and not is_biallelic_snp(rec):
            continue

        call = rec.samples[sample]
        code1 = get_phase_code(call)

        if code1 is None:
            continue

        n_vcf1_phased_het += 1

        ps = get_ps(call)

        if ps is None:
            n_no_ps += 1
            continue

        key = (rec.contig, rec.pos, rec.ref, rec.alts[0])

        if key not in vcf2_sites:
            continue

        code2 = vcf2_sites[key]["code"]

        # d = 0 means same orientation at this SNP.
        # d = 1 means flipped orientation at this SNP.
        d = code1 ^ code2

        ps_key = f"{rec.contig}:{ps}"

        phase_sets[ps_key].append(
            {
                "phase_set": ps_key,
                "chrom": rec.contig,
                "ps": ps,
                "pos": rec.pos,
                "ref": rec.ref,
                "alt": rec.alts[0],
                "vcf1_gt": "|".join(map(str, call.get("GT"))),
                "vcf2_gt": vcf2_sites[key]["gt"],
                "vcf1_code": code1,
                "vcf2_code": code2,
                "same_orientation": int(d == 0),
                "flipped_orientation": int(d == 1),
                "orientation_difference": d,
            }
        )

        n_shared += 1

    print(f"[VCF1] sample={sample}")
    print(f"[VCF1] records read: {n_total:,}")
    print(f"[VCF1] phased het SNPs: {n_vcf1_phased_het:,}")
    print(f"[VCF1] phased het SNPs without PS skipped: {n_no_ps:,}")
    print(f"[COMPARE] shared phased het SNPs with PS: {n_shared:,}")
    print(f"[COMPARE] phase sets with >=1 shared SNP: {len(phase_sets):,}")

    summary_rows = []
    snp_rows = []

    for ps_key, rows in phase_sets.items():
        rows = sorted(rows, key=lambda x: x["pos"])

        n = len(rows)
        n_same = sum(r["same_orientation"] for r in rows)
        n_flip = sum(r["flipped_orientation"] for r in rows)

        if n_same >= n_flip:
            best_orientation = "same"
            n_concordant_sites = n_same
            n_discordant_sites = n_flip
        else:
            best_orientation = "flipped"
            n_concordant_sites = n_flip
            n_discordant_sites = n_same

        site_concordance_after_best_flip = n_concordant_sites / n if n else math.nan

        # Pairwise phasing concordance:
        # pairs are concordant if both SNPs have the same orientation difference.
        total_pairs = n * (n - 1) // 2
        concordant_pairs = (n_same * (n_same - 1) // 2) + (n_flip * (n_flip - 1) // 2)

        if total_pairs > 0:
            pairwise_concordance = concordant_pairs / total_pairs
            pairwise_discordance = 1 - pairwise_concordance
        else:
            pairwise_concordance = math.nan
            pairwise_discordance = math.nan

        # Switch-like count across ordered SNPs:
        # how often the VCF2 orientation relative to VCF1 changes along the PS.
        diffs = [r["orientation_difference"] for r in rows]
        n_switches = sum(diffs[i] != diffs[i - 1] for i in range(1, len(diffs)))
        switch_rate = n_switches / (n - 1) if n > 1 else math.nan

        chrom = rows[0]["chrom"]
        ps = rows[0]["ps"]
        start = rows[0]["pos"]
        end = rows[-1]["pos"]

        summary_rows.append(
            {
                "phase_set": ps_key,
                "chrom": chrom,
                "ps": ps,
                "start": start,
                "end": end,
                "span_bp": end - start + 1,
                "n_shared_snps": n,
                "n_same_orientation": n_same,
                "n_flipped_orientation": n_flip,
                "best_orientation": best_orientation,
                "n_concordant_sites_after_best_flip": n_concordant_sites,
                "n_discordant_sites_after_best_flip": n_discordant_sites,
                "site_concordance_after_best_flip": site_concordance_after_best_flip,
                "total_snp_pairs": total_pairs,
                "concordant_snp_pairs": concordant_pairs,
                "pairwise_phase_concordance": pairwise_concordance,
                "pairwise_phase_discordance": pairwise_discordance,
                "n_orientation_switches": n_switches,
                "switch_rate": switch_rate,
            }
        )

        for r in rows:
            if best_orientation == "same":
                concordant_after_best_flip = r["orientation_difference"] == 0
            else:
                concordant_after_best_flip = r["orientation_difference"] == 1

            r["best_orientation_for_phase_set"] = best_orientation
            r["concordant_after_best_flip"] = int(concordant_after_best_flip)
            snp_rows.append(r)

    summary = pd.DataFrame(summary_rows)
    snps = pd.DataFrame(snp_rows)

    if not summary.empty:
        summary = summary.sort_values(
            ["chrom", "start", "end", "n_shared_snps"],
            ascending=[True, True, True, False],
        )

    if not snps.empty:
        snps = snps.sort_values(["chrom", "pos"])

    return summary, snps, sample


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare phasing concordance between two VCFs, using phase sets "
            "from VCF1 and allowing whole-phase-set flips."
        )
    )

    parser.add_argument("--vcf1", default=DEFAULT_VCF1, help="Reference phase-set VCF, e.g. WhatsHap VCF")
    parser.add_argument("--vcf2", default=DEFAULT_VCF2, help="Second phased VCF to compare against VCF1")
    parser.add_argument("--sample1", default=None, help="Sample name in VCF1. Default: first sample")
    parser.add_argument("--sample2", default=None, help="Sample name in VCF2. Default: first sample")
    parser.add_argument("--outdir", default="phase_concordance", help="Output directory")
    parser.add_argument(
        "--include-non-snps",
        action="store_true",
        help="Include non-SNP biallelic variants. Default: SNPs only.",
    )

    args = parser.parse_args()

    vcf1_path = Path(args.vcf1)
    vcf2_path = Path(args.vcf2)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    snps_only = not args.include_non_snps

    print("[INFO] VCF1:", vcf1_path)
    print("[INFO] VCF2:", vcf2_path)
    print("[INFO] SNPs only:", snps_only)

    vcf2_sites, sample2 = load_vcf2_sites(
        vcf2_path,
        sample=args.sample2,
        snps_only=snps_only,
    )

    summary, snps, sample1 = compare_phase_sets(
        vcf1_path,
        vcf2_sites,
        sample=args.sample1,
        snps_only=snps_only,
    )

    summary_file = outdir / "phase_set_concordance.tsv"
    snp_file = outdir / "per_snp_phase_comparison.tsv"

    summary.to_csv(summary_file, sep="\t", index=False)
    snps.to_csv(snp_file, sep="\t", index=False)

    print()
    print(f"[DONE] Wrote: {summary_file}")
    print(f"[DONE] Wrote: {snp_file}")

    if not summary.empty:
        print()
        print("[SUMMARY]")
        print(f"VCF1 sample: {sample1}")
        print(f"VCF2 sample: {sample2}")
        print(f"Phase sets compared: {len(summary):,}")
        print(f"Shared SNPs compared: {summary['n_shared_snps'].sum():,}")

        weighted_site_conc = (
            summary["site_concordance_after_best_flip"] * summary["n_shared_snps"]
        ).sum() / summary["n_shared_snps"].sum()

        pairwise_valid = summary.dropna(subset=["pairwise_phase_concordance"]).copy()

        if not pairwise_valid.empty and pairwise_valid["total_snp_pairs"].sum() > 0:
            weighted_pairwise_conc = (
                pairwise_valid["pairwise_phase_concordance"]
                * pairwise_valid["total_snp_pairs"]
            ).sum() / pairwise_valid["total_snp_pairs"].sum()
        else:
            weighted_pairwise_conc = math.nan

        print(f"Weighted site concordance after best PS flip: {weighted_site_conc:.4f}")
        print(f"Weighted pairwise phase concordance: {weighted_pairwise_conc:.4f}")

        print()
        print("[TOP 20 PHASE SETS BY N SHARED SNPS]")
        cols = [
            "phase_set",
            "start",
            "end",
            "n_shared_snps",
            "best_orientation",
            "site_concordance_after_best_flip",
            "pairwise_phase_concordance",
            "n_orientation_switches",
            "switch_rate",
        ]
        print(
            summary.sort_values("n_shared_snps", ascending=False)[cols]
            .head(20)
            .to_string(index=False)
        )


if __name__ == "__main__":
    main()