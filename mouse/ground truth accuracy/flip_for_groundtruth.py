#!/usr/bin/env python3

import gzip
import sys

in_vcf = sys.argv[1]
out_vcf = sys.argv[2]
def opener(path, mode):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

n_flipped = 0
n_total = 0

with opener(in_vcf, "rt") as fin, opener(out_vcf, "wt") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
            continue

        n_total += 1
        fields = line.rstrip("\n").split("\t")

        fmt = fields[8].split(":")
        samples = fields[9:]

        if "GT" not in fmt:
            fout.write(line)
            continue

        gt_i = fmt.index("GT")

        for i, sample in enumerate(samples):
            vals = sample.split(":")
            if len(vals) > gt_i and vals[gt_i] == "1|0":
                vals[gt_i] = "0|1"
                samples[i] = ":".join(vals)
                n_flipped += 1

        fields[9:] = samples
        fout.write("\t".join(fields) + "\n")

print(f"Variants processed: {n_total}")
print(f"Genotypes flipped 1|0 -> 0|1: {n_flipped}")
print(f"Wrote: {out_vcf}")