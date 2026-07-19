Simulated data uses the mouse per cell allele counts (available on figshare) as a starting point.

For all simulations, the following are the parameters of the simulations that were tested:
- escape: The percentage of genes that escape XCI. A gene specific Xi fraction is determined each gene selected to escape XCI, drawn randomly from a distribution from 0.1-0.5. Non escape genes have their Xi fraction drawn from a distribution from 0-0.1 Once the gene specific Xi fraction is selected, it is applied to all expressed SNPs of that gene.
- x1_fraction: The fraction of cells where X1 is the active X.
- snp_removal: The percentage of SNPs to remove from the per cell allele counts.
- cell_removal: The percentage of cells to remove from the per cell allele counts.
- read_removal: The percentage of allelic reads to remove from the per cell allele counts.

Each parameter was set as follows, unless it was the parameter of interest where a sweep from 0-100 was performed.
- escape: 25
- x1_fraction: 50
- snp_removal: 80
- cell_removal: 0
- read_removal: 0

Simulation commands were generated as follows e.g for cell removal:
```python
from pathlib import Path

out = Path("simulation_scDaisychain_master_cmds.txt")

escape = 25
snp_removal = 80
x1_fraction = 0.5
read_cutoff = 10
read_removal = 0

cell_removals = list(range(5, 96, 5)) + [99, 99.1, 99.2,99.3,99.4,99.5,99.6,99.7,99.8,99.9]
iterations = range(1, 201)

with out.open("w") as f:
    for cell_removal in cell_removals:
        for iteration in iterations:
            f.write(
                f"bash simulation_scDaisychain_master_run.sh "
                f"{escape:02d} {snp_removal} {iteration} "
                f"{x1_fraction:g} {read_cutoff} {cell_removal} {read_removal}\n"
            )

print(f"Wrote {out}")
print(f"Total commands: {len(cell_removals) * len(iterations)}")
print(f"Escape value: {escape}")
print(f"SNP removal value: {snp_removal}")
print(f"X1 fraction: {x1_fraction}")
print(f"Read cutoff: {read_cutoff}")
print(f"Read removal: {read_removal}")
print(f"Cell removal values: {cell_removals}")
```
A summary of the simulation output is provided on [figshare](https://figshare.com/account/articles/32942879). 

Accuracy was measured by the number of correctly phased SNPs / the total number of phased SNPs produced by running scDaisychain, using the workbook analysis.ipynb
