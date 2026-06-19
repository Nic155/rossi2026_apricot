#!/usr/bin/env python3
"""
Parse per-individual ROH intervals from a PLINK .hom file and intersect with CDS.

The .hom file is whitespace-delimited with a header. Relevant columns:
  IID  (column index 1) -- sample identifier
  CHR  (column index 3) -- chromosome as integer; converted to chr1, chr2 etc.
  POS1 (column index 6) -- ROH start, 1-based inclusive
  POS2 (column index 7) -- ROH end, 1-based inclusive

Coordinates are converted to 0-based half-open BED:
  BED start = POS1 - 1
  BED end   = POS2

For each individual, this script produces:
  results/regions/roh_per_sample/{sample}.bed      -- sorted ROH intervals
  results/regions/roh_per_sample_cds/{sample}.bed  -- ROH intersected with CDS
  results/regions/nonroh_per_sample_cds/{sample}.bed -- CDS minus ROH

A per-individual CDS-length summary is written to
  results/qc/per_individual_cds_lengths.tsv
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


WORKDIR = Path("~/work/ROH/roh_gl").expanduser()
DEFAULT_HOM = Path("~/work/ROH/roh/roh_800kb.hom").expanduser()
DEFAULT_CDS = Path("results/regions/cds_reduced.bed")


def fail(message):
    print("ERROR: {}".format(message), file=sys.stderr)
    sys.exit(1)


def require_file(path, label):
    if not path.exists() or path.stat().st_size == 0:
        fail("{} does not exist or is empty: {}".format(label, path))


def normalize_chr(value):
    chrom = str(value).strip()
    if not chrom or chrom.lower() == "nan":
        fail("Encountered an empty chromosome value in the .hom file.")
    if chrom.startswith("chr"):
        return chrom
    if chrom.endswith(".0"):
        chrom = chrom[:-2]
    return "chr{}".format(chrom)


def chrom_sort_key(chrom):
    raw = chrom[3:] if chrom.startswith("chr") else chrom
    try:
        return (0, int(raw))
    except ValueError:
        return (1, raw)


def bed_total_bp(path):
    total = 0
    if not path.exists():
        return 0
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            total += int(fields[2]) - int(fields[1])
    return total


def write_sorted_bed(intervals, path):
    intervals.sort(key=lambda item: (chrom_sort_key(item[0]), item[1], item[2]))
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(intervals)


def run_bedtools(args, output_path):
    with output_path.open("w") as out_handle:
        subprocess.run(args, check=True, stdout=out_handle)


def main():
    parser = argparse.ArgumentParser(
        description="Parse ROH intervals from .hom and intersect with CDS per individual."
    )
    parser.add_argument("--hom", default=str(DEFAULT_HOM), help="PLINK .hom file")
    parser.add_argument("--cds", default=str(DEFAULT_CDS), help="Reduced CDS BED")
    parser.add_argument(
        "--roh-dir",
        default="results/regions/roh_per_sample",
        help="Output directory for per-sample ROH BED files.",
    )
    parser.add_argument(
        "--roh-cds-dir",
        default="results/regions/roh_per_sample_cds",
        help="Output directory for per-sample ROH intersected with CDS.",
    )
    parser.add_argument(
        "--nonroh-cds-dir",
        default="results/regions/nonroh_per_sample_cds",
        help="Output directory for per-sample CDS minus ROH.",
    )
    parser.add_argument(
        "--qc",
        default="results/qc/per_individual_cds_lengths.tsv",
        help="Output QC table of per-individual CDS bp.",
    )
    args = parser.parse_args()

    os.chdir(str(WORKDIR))

    hom_path = Path(args.hom).expanduser()
    cds_path = Path(args.cds)
    roh_dir = Path(args.roh_dir)
    roh_cds_dir = Path(args.roh_cds_dir)
    nonroh_cds_dir = Path(args.nonroh_cds_dir)
    qc_path = Path(args.qc)

    require_file(hom_path, "PLINK .hom file")
    require_file(cds_path, "Reduced CDS BED")

    if shutil.which("bedtools") is None:
        fail("bedtools was not found in PATH. Load bedtools before running this step.")

    roh_dir.mkdir(parents=True, exist_ok=True)
    roh_cds_dir.mkdir(parents=True, exist_ok=True)
    nonroh_cds_dir.mkdir(parents=True, exist_ok=True)
    qc_path.parent.mkdir(parents=True, exist_ok=True)

    print("Reading ROH intervals from {}".format(hom_path))
    roh_by_sample = defaultdict(list)
    n_roh_total = 0

    with hom_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.split()
            if fields[0] == "FID":
                continue
            if len(fields) < 8:
                fail("Too few fields at line {}: {}".format(line_number, line.rstrip()))
            sample = fields[1]
            try:
                chrom = normalize_chr(fields[3])
                start = int(fields[6]) - 1
                end = int(fields[7])
            except (ValueError, IndexError):
                fail("Could not parse ROH interval at line {}: {}".format(line_number, line.rstrip()))
            if start < 0 or end <= start:
                fail(
                    "Invalid BED coordinates at line {}: start={} end={}".format(
                        line_number, start, end
                    )
                )
            roh_by_sample[sample].append((chrom, start, end))
            n_roh_total += 1

    n_samples = len(roh_by_sample)
    print("Samples with ROH intervals: {}".format(n_samples))
    print("Total ROH intervals: {}".format(n_roh_total))

    qc_rows = []

    for sample, intervals in sorted(roh_by_sample.items()):
        roh_bed = roh_dir / "{}.bed".format(sample)
        roh_cds_bed = roh_cds_dir / "{}.bed".format(sample)
        nonroh_cds_bed = nonroh_cds_dir / "{}.bed".format(sample)

        write_sorted_bed(intervals, roh_bed)

        run_bedtools(
            ["bedtools", "intersect", "-a", str(cds_path), "-b", str(roh_bed)],
            roh_cds_bed,
        )
        run_bedtools(
            ["bedtools", "subtract", "-a", str(cds_path), "-b", str(roh_bed)],
            nonroh_cds_bed,
        )

        roh_cds_bp = bed_total_bp(roh_cds_bed)
        nonroh_cds_bp = bed_total_bp(nonroh_cds_bed)
        qc_rows.append((sample, len(intervals), roh_cds_bp, nonroh_cds_bp))

    print("Writing per-individual CDS lengths to {}".format(qc_path))
    with qc_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["genotype", "n_roh_intervals", "roh_cds_bp", "nonroh_cds_bp"])
        for row in qc_rows:
            writer.writerow(row)

    total_roh_cds_bp = sum(r[2] for r in qc_rows)
    print("Samples processed: {}".format(len(qc_rows)))
    print("Total ROH CDS bp across all samples: {}".format(total_roh_cds_bp))
    print("Wrote per-sample ROH BEDs to {}".format(roh_dir))
    print("Wrote per-sample ROH CDS BEDs to {}".format(roh_cds_dir))
    print("Wrote per-sample non-ROH CDS BEDs to {}".format(nonroh_cds_dir))
    print("Wrote {}".format(qc_path))


if __name__ == "__main__":
    main()
