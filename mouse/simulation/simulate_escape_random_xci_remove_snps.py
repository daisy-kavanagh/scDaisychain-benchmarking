#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument(
        "--counts",
        default="/home/913/dk4874/scratch/gdata/scDaisychain_paper/mouse/output/variant_counts_fast_filtered_readids_from_WGS/mouse_per_cell_allele_counts.filtered.tsv",
    )
    p.add_argument(
        "--active-x",
        default="/home/913/dk4874/scratch/gdata/scDaisychain_paper/mouse/output/scDaisychain_modes_multi_filtered_weighted_from_WGS/groundtruth_split/haplotype_sums_on_the_fly.csv",
    )
    p.add_argument(
        "--outdir",
        default="/home/913/dk4874/scratch/gdata/scDaisychain_paper/mouse/scripts/simulation/simulated_escape_counts_from_WGS",
    )

    p.add_argument("--escape-pct", type=int, required=True)
    p.add_argument("--snp-removal-pct", type=float, required=True)
    p.add_argument("--iteration", type=int, required=True)

    p.add_argument("--seed", type=int, default=1)

    p.add_argument("--x1-fraction", type=float, default=0.5)
    p.add_argument("--cell-removal-pct", type=float, default=0.0)
    p.add_argument("--read-removal-pct", type=float, default=0.0)

    p.add_argument("--escape-min", type=float, default=0.10)
    p.add_argument("--escape-max", type=float, default=0.50)
    p.add_argument("--nonescape-min", type=float, default=0.00)
    p.add_argument("--nonescape-max", type=float, default=0.00)

    return p.parse_args()


def clean_gene(x):
    if pd.isna(x):
        return ""
    s = str(x).strip()
    return "" if s in {"", "-", ".", "NA", "NAN", "nan", "None"} else s


def clean_cb(x):
    s = str(x).strip()
    return s[:-2] if s.endswith("-1") else s


def label_float(prefix, x):
    return f"{prefix}{x:g}".replace(".", "p")


def main():
    args = parse_args()

    if not 0 <= args.x1_fraction <= 1:
        raise ValueError("--x1-fraction must be between 0 and 1")

    if not 0 <= args.cell_removal_pct <= 100:
        raise ValueError("--cell-removal-pct must be between 0 and 100")

    if not 0 <= args.read_removal_pct <= 100:
        raise ValueError("--read-removal-pct must be between 0 and 100")

    seed = (
        int(args.seed)
        + int(args.iteration) * 1_000_000
        + int(args.escape_pct) * 10_000
        + int(round(args.snp_removal_pct * 100))
        + int(round(args.x1_fraction * 1000))
        + int(round(args.cell_removal_pct * 100)) * 10
        + int(round(args.read_removal_pct * 100))
    )

    rng = np.random.default_rng(seed)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    escape_label = f"escape{args.escape_pct:02d}"
    snp_label = label_float("snpremove", args.snp_removal_pct)
    skew_label = label_float("skewX1", args.x1_fraction)
    cell_label = label_float("cellremove", args.cell_removal_pct)
    readremove_label = label_float("readremove", args.read_removal_pct)
    iter_label = f"iter{args.iteration:03d}"

    prefix = (
        f"{escape_label}.{snp_label}.{skew_label}."
        f"{cell_label}.{readremove_label}.{iter_label}"
    )

    x = pd.read_csv(args.counts, sep="\t")
    ax = pd.read_csv(args.active_x)

    required = {"cell_barcode", "active_X"}
    if not required.issubset(ax.columns):
        raise ValueError("--active-x must contain columns: cell_barcode, active_X")

    for col in ["cell_barcode", "gene", "refCount", "altCount", "otherCount"]:
        if col not in x.columns:
            raise ValueError(f"--counts is missing required column: {col}")

    if "position" not in x.columns:
        raise ValueError("--counts must contain a SNP position column called 'position'")

    if "contig" not in x.columns:
        x["contig"] = "chrX"

    x["cell_barcode"] = x["cell_barcode"].astype(str).map(clean_cb)
    ax["cell_barcode"] = ax["cell_barcode"].astype(str).map(clean_cb)

    active_map = ax.drop_duplicates("cell_barcode").set_index("cell_barcode")["active_X"]

    x["original_active_X"] = x["cell_barcode"].map(active_map)
    missing = x["original_active_X"].isna().sum()
    if missing:
        print(f"[WARN] rows with no original active_X call: {missing:,}; dropping")
        x = x.dropna(subset=["original_active_X"]).copy()

    x = x[x["original_active_X"].isin(["X1", "X2"])].copy()
    x["gene_clean"] = x["gene"].map(clean_gene)

    # Remove cells before SNP/read simulation.
    cells_initial = np.array(sorted(x["cell_barcode"].unique()))
    n_cells_initial = len(cells_initial)

    n_cells_remove = int(round(n_cells_initial * args.cell_removal_pct / 100))

    if n_cells_remove > 0:
        removed_cells = set(
            rng.choice(cells_initial, size=n_cells_remove, replace=False)
        )
    else:
        removed_cells = set()

    x["cell_removed"] = x["cell_barcode"].isin(removed_cells)
    x = x.loc[~x["cell_removed"]].copy()

    cells_after_cell_removal = np.array(sorted(x["cell_barcode"].unique()))
    n_cells_after_cell_removal = len(cells_after_cell_removal)

    removed_cells_df = pd.DataFrame({
        "cell_barcode": sorted(removed_cells),
        "removed": True,
        "cell_removal_pct_setting": args.cell_removal_pct,
        "escape_pct_setting": args.escape_pct,
        "snp_removal_pct_setting": args.snp_removal_pct,
        "x1_fraction_setting": args.x1_fraction,
        "read_removal_pct_setting": args.read_removal_pct,
        "iteration": args.iteration,
    })

    # Define SNPs by genomic position.
    snp_cols = ["contig", "position"]
    snps = x[snp_cols].drop_duplicates().copy()
    n_snps_initial = len(snps)

    n_remove = int(round(n_snps_initial * args.snp_removal_pct / 100))

    if n_remove > 0:
        remove_idx = rng.choice(snps.index.to_numpy(), size=n_remove, replace=False)
        removed_snps = snps.loc[remove_idx].copy()
    else:
        removed_snps = snps.iloc[0:0].copy()

    removed_snps["removed"] = True

    x = x.merge(
        removed_snps[snp_cols + ["removed"]],
        on=snp_cols,
        how="left",
    )
    x["removed"] = x["removed"].fillna(False)

    y = x.loc[~x["removed"]].copy()

    n_snps_after = y[snp_cols].drop_duplicates().shape[0]
    genes_initial = sorted(x.loc[x["gene_clean"] != "", "gene_clean"].unique())
    genes_after = sorted(y.loc[y["gene_clean"] != "", "gene_clean"].unique())

    n_escape = int(round(len(genes_after) * args.escape_pct / 100))
    escape_genes = (
        set(rng.choice(genes_after, size=n_escape, replace=False))
        if n_escape
        else set()
    )

    xi_frac_by_gene = {}
    gene_meta = []

    for g in genes_after:
        is_escape = g in escape_genes

        if is_escape:
            xi_frac = rng.uniform(args.escape_min, args.escape_max)
        else:
            xi_frac = rng.uniform(args.nonescape_min, args.nonescape_max)

        xi_frac_by_gene[g] = xi_frac

        gene_rows = y["gene_clean"].eq(g)

        gene_meta.append({
            "gene": g,
            "is_escape_gene": is_escape,
            "target_xi_fraction": xi_frac,
            "escape_pct_setting": args.escape_pct,
            "snp_removal_pct_setting": args.snp_removal_pct,
            "x1_fraction_setting": args.x1_fraction,
            "cell_removal_pct_setting": args.cell_removal_pct,
            "read_removal_pct_setting": args.read_removal_pct,
            "iteration": args.iteration,
            "n_rows_after_snp_removal": int(gene_rows.sum()),
            "n_snps_after_snp_removal": int(
                y.loc[gene_rows, snp_cols].drop_duplicates().shape[0]
            ),
        })

    allele_total_before_read_removal = (
        y["refCount"].astype(int)
        + y["altCount"].astype(int)
    )

    keep_read_prob = 1.0 - (args.read_removal_pct / 100.0)

    allele_total = rng.binomial(
        allele_total_before_read_removal.to_numpy(),
        keep_read_prob,
    )

    xi_frac = y["gene_clean"].map(xi_frac_by_gene).fillna(0.0).astype(float)

    xi_reads = rng.binomial(
        allele_total,
        xi_frac.to_numpy(),
    )
    xa_reads = allele_total - xi_reads

    # Simulated XCI skew.
    cells = np.array(sorted(y["cell_barcode"].unique()))

    sim_active_x = rng.choice(
        ["X1", "X2"],
        size=len(cells),
        replace=True,
        p=[args.x1_fraction, 1 - args.x1_fraction],
    )

    cell_active_df = pd.DataFrame({
        "cell_barcode": cells,
        "sim_active_X": sim_active_x,
        "x1_fraction_setting": args.x1_fraction,
        "cell_removal_pct_setting": args.cell_removal_pct,
        "read_removal_pct_setting": args.read_removal_pct,
        "escape_pct_setting": args.escape_pct,
        "snp_removal_pct_setting": args.snp_removal_pct,
        "iteration": args.iteration,
    })

    sim_active_map = cell_active_df.set_index("cell_barcode")["sim_active_X"]
    y["sim_active_X"] = y["cell_barcode"].map(sim_active_map)

    is_x1 = y["sim_active_X"].eq("X1").to_numpy()

    # Convention:
    # H1/X1 = ALT
    # H2/X2 = REF
    #
    # X1 active: ALT = Xa, REF = Xi
    # X2 active: REF = Xa, ALT = Xi
    y["refCount"] = np.where(is_x1, xi_reads, xa_reads).astype(int)
    y["altCount"] = np.where(is_x1, xa_reads, xi_reads).astype(int)
    y["totalCount"] = (
        y["refCount"]
        + y["altCount"]
        + y["otherCount"].astype(int)
    )

    removed_snps.to_csv(
        outdir / f"removed_snps.{prefix}.tsv",
        sep="\t",
        index=False,
    )

    removed_cells_df.to_csv(
        outdir / f"removed_cells.{prefix}.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(gene_meta).to_csv(
        outdir / f"simulation_gene_metadata.{prefix}.tsv",
        sep="\t",
        index=False,
    )

    cell_active_df.to_csv(
        outdir / f"simulation_cell_active_x.{prefix}.tsv",
        sep="\t",
        index=False,
    )

    cell_summary = (
        cell_active_df["sim_active_X"]
        .value_counts()
        .rename_axis("sim_active_X")
        .reset_index(name="n_cells")
    )

    cell_summary["fraction_cells"] = (
        cell_summary["n_cells"] / cell_summary["n_cells"].sum()
    )

    cell_summary["x1_fraction_setting"] = args.x1_fraction
    cell_summary["cell_removal_pct_setting"] = args.cell_removal_pct
    cell_summary["read_removal_pct_setting"] = args.read_removal_pct
    cell_summary["escape_pct_setting"] = args.escape_pct
    cell_summary["snp_removal_pct_setting"] = args.snp_removal_pct
    cell_summary["iteration"] = args.iteration

    cell_summary.to_csv(
        outdir / f"simulation_cell_active_x_summary.{prefix}.tsv",
        sep="\t",
        index=False,
    )

    summary = pd.DataFrame([{
        "escape_pct_setting": args.escape_pct,
        "snp_removal_pct_setting": args.snp_removal_pct,
        "x1_fraction_setting": args.x1_fraction,
        "cell_removal_pct_setting": args.cell_removal_pct,
        "read_removal_pct_setting": args.read_removal_pct,
        "iteration": args.iteration,
        "seed": seed,
        "n_rows_initial": len(x),
        "n_rows_after_snp_removal": len(y),
        "n_cells_initial": n_cells_initial,
        "n_cells_removed": n_cells_remove,
        "n_cells_after_cell_removal": n_cells_after_cell_removal,
        "n_snps_initial": n_snps_initial,
        "n_snps_removed": n_remove,
        "n_snps_after_snp_removal": n_snps_after,
        "allele_reads_before_read_removal": int(allele_total_before_read_removal.sum()),
        "allele_reads_after_read_removal": int(allele_total.sum()),
        "n_genes_initial": len(genes_initial),
        "n_genes_after_snp_removal": len(genes_after),
        "n_escape_genes": n_escape,
        "n_nonescape_genes": len(genes_after) - n_escape,
        "escape_min": args.escape_min,
        "escape_max": args.escape_max,
        "nonescape_min": args.nonescape_min,
        "nonescape_max": args.nonescape_max,
    }])

    summary.to_csv(
        outdir / f"simulation_summary.{prefix}.tsv",
        sep="\t",
        index=False,
    )

    drop_cols = [
        "gene_clean",
        "original_active_X",
        "removed",
        "cell_removed",
        "sim_active_X",
    ]

    y = y.drop(columns=[c for c in drop_cols if c in y.columns])

    out_counts = outdir / f"mouse_per_cell_allele_counts.{prefix}.tsv"
    y.to_csv(out_counts, sep="\t", index=False)

    print(f"[done] wrote counts: {out_counts}")
    print(f"[done] wrote summary: {outdir / f'simulation_summary.{prefix}.tsv'}")
    print(
        f"[done] simulated X1 fraction setting: {args.x1_fraction}; "
        f"observed: {cell_active_df['sim_active_X'].eq('X1').mean():.4f}"
    )
    print(
        f"[done] cell removal setting: {args.cell_removal_pct}%; "
        f"removed {n_cells_remove:,}/{n_cells_initial:,} cells"
    )
    print(
        f"[done] read removal setting: {args.read_removal_pct}%; "
        f"kept {int(allele_total.sum()):,}/"
        f"{int(allele_total_before_read_removal.sum()):,} allele reads"
    )


if __name__ == "__main__":
    main()