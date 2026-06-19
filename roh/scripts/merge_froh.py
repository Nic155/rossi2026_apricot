#!/usr/bin/env python3
"""Parse per-threshold .hom.indiv files and compute F_ROH for each sample."""

import os
import sys

WORKDIR = "/home/nrossi/work/ROH"
GENOME_SIZE_KB = 220000.0  # 220 Mb expressed in kb

THRESHOLDS = [100, 500, 800, 1000, 2000]


def parse_hom_indiv(path):
    """Return dict {IID: total_KB} from a .hom.indiv file."""
    result = {}
    with open(path) as fh:
        header = fh.readline().split()
        iid_idx = header.index("IID")
        kb_idx = header.index("KB")
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            result[parts[iid_idx]] = float(parts[kb_idx])
    return result


data = {}  # {sample: {threshold_kb: froh}}

for t in THRESHOLDS:
    fname = os.path.join(WORKDIR, "roh", f"roh_{t}kb.hom.indiv")
    if not os.path.isfile(fname):
        sys.exit(f"ERROR: expected file not found: {fname}")
    per_sample = parse_hom_indiv(fname)
    for iid, kb in per_sample.items():
        data.setdefault(iid, {})[t] = kb / GENOME_SIZE_KB

outdir = os.path.join(WORKDIR, "results")
os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, "FROH_all_thresholds.tsv")

col_names = [f"froh_{t}kb" for t in THRESHOLDS]
with open(outfile, "w") as fh:
    fh.write("sample\t" + "\t".join(col_names) + "\n")
    for sample in sorted(data.keys()):
        vals = []
        for t in THRESHOLDS:
            v = data[sample].get(t)
            vals.append(f"{v:.6f}" if v is not None else "NA")
        fh.write(sample + "\t" + "\t".join(vals) + "\n")

print(f"Written {len(data)} samples to {outfile}")
