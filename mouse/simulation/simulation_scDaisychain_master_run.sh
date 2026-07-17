#!/bin/bash

escape_percentage=$1
snp_removal_percentage=$2
iteration=$3

# optional args
x1_fraction=${4:-0.5}
read_cutoff=${5:-10}
cell_removal_percentage=${6:-0}
read_removal_percentage=${7:-0}

escape_label=$(printf "escape%02d" "$escape_percentage")

snp_label=$(python - <<PY
x = float("$snp_removal_percentage")
print(f"snpremove{x:g}".replace(".", "p"))
PY
)

iter_label=$(printf "iter%03d" "$iteration")

skew_label=$(python - <<PY
x = float("$x1_fraction")
print(f"skewX1{x:g}".replace(".", "p"))
PY
)

cell_label=$(python - <<PY
x = float("$cell_removal_percentage")
print(f"cellremove{x:g}".replace(".", "p"))
PY
)

readremove_label=$(python - <<PY
x = float("$read_removal_percentage")
print(f"readremove{x:g}".replace(".", "p"))
PY
)

read_label=$(printf "minreads%02d" "$read_cutoff")

base_sim="/home/913/dk4874/scratch/vn68_gdata/scDaisychain_simulations/snp_removal2"
sim_counts_dir="${base_sim}/simulated_escape_counts_from_WGS"

outdir="${base_sim}/scDaisychain/${escape_label}/${snp_label}/${skew_label}/${cell_label}/${readremove_label}/${read_label}/${iter_label}"
mkdir -p "$outdir"

echo "Running:"
echo "  escape_percentage=${escape_percentage}"
echo "  snp_removal_percentage=${snp_removal_percentage}"
echo "  iteration=${iteration}"
echo "  x1_fraction=${x1_fraction}"
echo "  read_cutoff=${read_cutoff}"
echo "  cell_removal_percentage=${cell_removal_percentage}"
echo "  read_removal_percentage=${read_removal_percentage}"
echo "  outdir=${outdir}"

source ~/.bashrc
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

module load htslib || true

###############################################################################
# 1. Simulate counts
###############################################################################

mamba activate X_Inactivation2

python simulate_escape_random_xci_remove_snps.py \
  --escape-pct "$escape_percentage" \
  --snp-removal-pct "$snp_removal_percentage" \
  --iteration "$iteration" \
  --x1-fraction "$x1_fraction" \
  --cell-removal-pct "$cell_removal_percentage" \
  --read-removal-pct "$read_removal_percentage" \
  --outdir "$sim_counts_dir"

counts="${sim_counts_dir}/mouse_per_cell_allele_counts.${escape_label}.${snp_label}.${skew_label}.${cell_label}.${readremove_label}.${iter_label}.tsv"

###############################################################################
# 2. Run scDaisychain phasing
###############################################################################

python scDaisychain-phase-x \
  "$counts" \
  "$outdir" \
  --min_reads "$read_cutoff" \
  --lower_cutoff 0.01 \
  --partition_mode weighted \
  --stage2-mode cell