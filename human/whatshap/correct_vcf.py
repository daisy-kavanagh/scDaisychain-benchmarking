#!/usr/bin/env python3

import argparse
import gzip
import pandas as pd
import pysam


parser = argparse.ArgumentParser(
    description=(
        "Flip discordant SNP genotypes in a DaisyChain VCF according to "
        "per_snp_phase_comparison.tsv majority phase-set orientation."
    )
)

parser.add_argument("--daisychain-vcf", required=True)
parser.add_argument("--comparison-tsv", required=True)
parser.add_argument("--out-vcf", required=True)

parser.add_argument(
    "--chrom-col",
    default="chrom",
    help="Chromosome column in comparison TSV"
)

parser.add_argument(
    "--pos-col",
    default="pos",
    help="Position column in comparison TSV"
)

parser.add_argument(
    "--best-col",
    default="best_orientation_for_phase_set",
    help="Column giving phase-set best orientation: same/flipped"
)

parser.add_argument(
    "--same-col",
    default="same_orientation",
    help="Column indicating SNP is same orientation"
)

parser.add_argument(
    "--flipped-col",
    default="flipped_orientation",
    help="Column indicating SNP is flipped orientation"
)

parser.add_argument(
    "--only-discordant-after-best-flip",
    action="store_true",
    help=(
        "If set, only consider rows where concordant_after_best_flip == 0. "
        "Normally this is equivalent to the same/flipped logic."
    )
)

args = parser.parse_args()


def as_int01(x):
    if pd.isna(x):
        return 0
    if isinstance(x, bool):
        return int(x)
    s = str(x).strip().lower()
    return 1 if s in {"1", "true", "t", "yes", "y"} else 0


def flip_gt_tuple(gt):
    """
    Flip phased genotype orientation.
    0|1 -> 1|0
    1|0 -> 0|1
    0|0 -> 0|0
    1|1 -> 1|1
    """
    if gt is None:
        return gt
    if len(gt) != 2:
        return gt
    if gt[0] is None or gt[1] is None:
        return gt
    return (gt[1], gt[0])


# ------------------------------------------------------------
# Read comparison table and identify SNPs to flip in DaisyChain
# ------------------------------------------------------------

comp = pd.read_csv(args.comparison_tsv, sep="\t")

required = {
    args.chrom_col,
    args.pos_col,
    args.best_col,
    args.same_col,
    args.flipped_col,
}

missing = required - set(comp.columns)
if missing:
    raise ValueError(f"Missing required columns in comparison TSV: {missing}")

if args.only_discordant_after_best_flip:
    if "concordant_after_best_flip" not in comp.columns:
        raise ValueError(
            "--only-discordant-after-best-flip requires concordant_after_best_flip column"
        )
    comp = comp[comp["concordant_after_best_flip"].map(as_int01).eq(0)].copy()

to_flip = set()

for _, row in comp.iterrows():
    chrom = str(row[args.chrom_col])
    pos = int(row[args.pos_col])

    best = str(row[args.best_col]).strip().lower()
    same = as_int01(row[args.same_col])
    flipped = as_int01(row[args.flipped_col])

    flip_this = False

    if best == "same":
        # Majority of phase set is same, so flipped SNPs are discordant.
        if flipped == 1:
            flip_this = True

    elif best == "flipped":
        # Majority of phase set is flipped, so same-orientation SNPs are discordant.
        if same == 1:
            flip_this = True

    else:
        continue

    if flip_this:
        to_flip.add((chrom, pos))

print(f"[INFO] SNPs marked for DaisyChain GT flipping: {len(to_flip):,}")

# ------------------------------------------------------------
# Write list of SNPs selected for flipping
# ------------------------------------------------------------

flip_sites_tsv = args.out_vcf + ".flipped_sites.tsv"

with open(flip_sites_tsv, "w") as f:
    f.write("chrom\tpos\n")
    for chrom, pos in sorted(to_flip):
        f.write(f"{chrom}\t{pos}\n")

print(f"[INFO] Wrote flipped-site list: {flip_sites_tsv}")


# ------------------------------------------------------------
# Edit VCF
# ------------------------------------------------------------

vcf_in = pysam.VariantFile(args.daisychain_vcf)
vcf_out = pysam.VariantFile(args.out_vcf, "wz", header=vcf_in.header)

n_records = 0
n_flipped = 0
n_target_sites_seen = 0
n_skipped_unphased = 0
n_skipped_missing_gt = 0

samples = list(vcf_in.header.samples)

if len(samples) != 1:
    print(
        f"[WARN] VCF has {len(samples)} samples. "
        "This script will flip GT for all samples."
    )

for rec in vcf_in.fetch():
    n_records += 1

    key = (rec.chrom, rec.pos)

    if key in to_flip:
        n_target_sites_seen += 1
        did_flip = False

        for sample in samples:
            sample_data = rec.samples[sample]
            gt = sample_data.get("GT")

            if gt is None or len(gt) != 2:
                n_skipped_missing_gt += 1
                continue

            if not sample_data.phased:
                n_skipped_unphased += 1
                continue

            # Flip phase orientation:
            # 0|1 -> 1|0
            # 1|0 -> 0|1
            sample_data["GT"] = (gt[1], gt[0])
            sample_data.phased = True
            did_flip = True

        if did_flip:
            n_flipped += 1

    vcf_out.write(rec)

vcf_in.close()
vcf_out.close()

print(f"[INFO] VCF records processed: {n_records:,}")
print(f"[INFO] Target flip sites from TSV: {len(to_flip):,}")
print(f"[INFO] Target flip sites seen in VCF: {n_target_sites_seen:,}")
print(f"[INFO] VCF records actually flipped: {n_flipped:,}")
print(f"[INFO] Missing-GT target records skipped: {n_skipped_missing_gt:,}")
print(f"[INFO] Unphased target records skipped: {n_skipped_unphased:,}")
print(f"[INFO] Wrote corrected VCF: {args.out_vcf}")