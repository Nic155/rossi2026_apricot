#!/usr/bin/env python3
"""
Calculate per-genotype deleterious dosage per 100 kb CDS inside and outside ROH.

For each individual the pipeline uses per-individual ROH CDS intervals produced
by step 02. The load variant BED from step 03 is intersected with each
individual's ROH regions using bedtools to classify variants as inside or
outside ROH CDS.

Formula:
  roh_deleterious_dosage_per_100kb_CDS =
      sum dosage for PROVEAN/HIGH variants in that individual's ROH CDS
      / that individual's ROH CDS bp * 100000

  nonroh_deleterious_dosage_per_100kb_CDS =
      sum dosage for PROVEAN/HIGH variants in that individual's non-ROH CDS
      / that individual's non-ROH CDS bp * 100000

NA genotype calls are ignored in dosage sums. NA is not treated as zero.

Individuals with no ROH intervals (absent from the .hom file or with empty ROH
BED) have roh_deleterious_dosage_per_100kb_CDS reported as NA. Their
nonroh_cds_bp equals the total reduced CDS bp.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path


WORKDIR = Path("~/work/ROH/roh_gl").expanduser()
DEFAULT_POPMAP = Path("~/work/ROH/list.txt").expanduser()
META_COLS = ["CHR", "POS", "AA", "IMPACT", "EFFECT", "PROVEAN_SCORE", "DEL_CLASS"]
MISSING_VALUES = {"", "NA", "NAN", "NULL", "."}


def fail(message):
    print("ERROR: {}".format(message), file=sys.stderr)
    sys.exit(1)


def require_file(path, label):
    if not path.exists() or path.stat().st_size == 0:
        fail("{} does not exist or is empty: {}".format(label, path))


def read_popmap(path):
    popmap = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail("Population map has no header.")
        if "genotype" not in reader.fieldnames or "group" not in reader.fieldnames:
            fail("Population map must have 'genotype' and 'group' columns.")
        for row in reader:
            genotype = row.get("genotype", "").strip()
            group = row.get("group", "").strip()
            if genotype:
                popmap[genotype] = group
    return popmap


def read_total_cds_bp(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail("CDS lengths file has no header.")
        for row in reader:
            if row.get("region_class") == "all_CDS":
                return int(float(row["cds_length_bp"]))
    fail("CDS lengths file is missing region_class 'all_CDS'.")


def read_per_individual_lengths(path):
    lengths = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail("Per-individual CDS lengths file has no header.")
        required = {"genotype", "roh_cds_bp", "nonroh_cds_bp"}
        missing = required.difference(set(reader.fieldnames))
        if missing:
            fail(
                "Per-individual CDS lengths file is missing column(s): {}".format(
                    ", ".join(sorted(missing))
                )
            )
        for row in reader:
            genotype = row["genotype"]
            roh_cds_bp = int(row["roh_cds_bp"])
            nonroh_cds_bp = int(row["nonroh_cds_bp"])
            lengths[genotype] = (roh_cds_bp, nonroh_cds_bp)
    return lengths


def get_genotype_columns(path):
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            fail("{} is empty.".format(path))
    missing = [col for col in META_COLS if col not in header]
    if missing:
        fail("{} is missing metadata column(s): {}".format(path, ", ".join(missing)))
    genotype_cols = [col for col in header if col not in META_COLS]
    if not genotype_cols:
        fail("No genotype columns found in {}.".format(path))
    return genotype_cols


def parse_dosage(value, path, line_number, genotype):
    text = str(value).strip()
    if text.upper() in MISSING_VALUES:
        return None
    try:
        return float(text)
    except ValueError:
        fail(
            "Non-numeric non-missing genotype value '{}' for {} at {} line {}".format(
                value, genotype, path, line_number
            )
        )


def load_variant_dosages(tsv_path, bed_path, genotype_cols):
    """Return a dict mapping variant_id to a per-genotype dosage dict.

    Variant IDs are read from bed_path (col 4) and matched to TSV rows by
    position in file order, which is guaranteed to be identical by step 03.
    """
    variant_ids = []
    with bed_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                fail("Malformed BED line in {}: {}".format(bed_path, line.rstrip()))
            variant_ids.append(fields[3])

    variants = {}
    with tsv_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail("{} has no header.".format(tsv_path))
        missing = [col for col in genotype_cols if col not in reader.fieldnames]
        if missing:
            fail(
                "{} is missing genotype column(s): {}".format(
                    tsv_path, ", ".join(missing)
                )
            )
        for line_number, row in enumerate(reader, start=2):
            idx = line_number - 2
            if idx >= len(variant_ids):
                fail(
                    "TSV has more rows than BED has variant IDs. "
                    "Re-run step 03 to regenerate load_variants_cds files."
                )
            variant_id = variant_ids[idx]
            dosages = {}
            for genotype in genotype_cols:
                dosages[genotype] = parse_dosage(
                    row.get(genotype, ""), tsv_path, line_number, genotype
                )
            variants[variant_id] = dosages

    if len(variants) != len(variant_ids):
        fail(
            "TSV has fewer rows ({}) than BED has variant IDs ({}). "
            "Re-run step 03.".format(len(variants), len(variant_ids))
        )
    return variants


def get_roh_variant_ids(variant_bed, roh_bed):
    """Return the set of variant_ids from variant_bed that overlap roh_bed."""
    if not roh_bed.exists() or roh_bed.stat().st_size == 0:
        return set()
    result = subprocess.run(
        ["bedtools", "intersect", "-u", "-a", str(variant_bed), "-b", str(roh_bed)],
        check=True,
        capture_output=True,
        text=True,
    )
    ids = set()
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            fail("Malformed bedtools output line: {}".format(line.rstrip()))
        ids.add(fields[3])
    return ids


def main():
    parser = argparse.ArgumentParser(
        description="Calculate per-genotype ROH vs non-ROH deleterious burden."
    )
    parser.add_argument(
        "--variants-tsv", default="results/variants/load_variants_cds.tsv"
    )
    parser.add_argument(
        "--variants-bed", default="results/variants/load_variants_cds.bed"
    )
    parser.add_argument(
        "--cds-lengths", default="results/qc/cds_lengths.tsv"
    )
    parser.add_argument(
        "--per-individual-lengths",
        default="results/qc/per_individual_cds_lengths.tsv",
    )
    parser.add_argument(
        "--roh-dir", default="results/regions/roh_per_sample"
    )
    parser.add_argument(
        "--popmap", default=str(DEFAULT_POPMAP), help="Sample-population map TSV"
    )
    parser.add_argument(
        "--final-output",
        default="results/final_genotype_roh_nonroh_burden.tsv",
    )
    parser.add_argument(
        "--extended-qc", default="results/qc/genotype_burden_extended.tsv"
    )
    args = parser.parse_args()

    os.chdir(str(WORKDIR))

    variants_tsv = Path(args.variants_tsv)
    variants_bed = Path(args.variants_bed)
    cds_lengths = Path(args.cds_lengths)
    per_individual_lengths = Path(args.per_individual_lengths)
    roh_dir = Path(args.roh_dir)
    popmap_path = Path(args.popmap).expanduser()
    final_output = Path(args.final_output)
    extended_qc = Path(args.extended_qc)

    require_file(variants_tsv, "Load variants CDS TSV")
    require_file(variants_bed, "Load variants CDS BED")
    require_file(cds_lengths, "CDS lengths table")
    require_file(per_individual_lengths, "Per-individual CDS lengths table")
    require_file(popmap_path, "Sample-population map")

    if shutil.which("bedtools") is None:
        fail("bedtools was not found in PATH. Load bedtools before running this step.")

    final_output.parent.mkdir(parents=True, exist_ok=True)
    extended_qc.parent.mkdir(parents=True, exist_ok=True)

    total_cds_bp = read_total_cds_bp(cds_lengths)
    individual_lengths = read_per_individual_lengths(per_individual_lengths)
    popmap = read_popmap(popmap_path)
    genotype_cols = get_genotype_columns(variants_tsv)

    print("Genotype count: {}".format(len(genotype_cols)))
    print("Loading variant dosages from {}".format(variants_tsv))
    variants = load_variant_dosages(variants_tsv, variants_bed, genotype_cols)
    n_variants = len(variants)
    print("CDS load variants: {}".format(n_variants))
    if n_variants == 0:
        fail("No CDS load variants are present.")

    all_variant_ids = list(variants.keys())

    final_header = [
        "genotype",
        "population",
        "roh_deleterious_dosage_per_100kb_CDS",
        "nonroh_deleterious_dosage_per_100kb_CDS",
        "roh_cds_bp",
        "nonroh_cds_bp",
    ]
    extended_header = [
        "genotype",
        "population",
        "roh_raw_dosage_sum",
        "nonroh_raw_dosage_sum",
        "roh_cds_bp",
        "nonroh_cds_bp",
        "roh_n_variants",
        "nonroh_n_variants",
        "n_cds_variants_total",
        "roh_n_nonmissing_calls",
        "nonroh_n_nonmissing_calls",
        "roh_deleterious_dosage_per_100kb_CDS",
        "nonroh_deleterious_dosage_per_100kb_CDS",
    ]

    with final_output.open("w", newline="") as final_handle, \
            extended_qc.open("w", newline="") as extended_handle:
        final_writer = csv.writer(final_handle, delimiter="\t", lineterminator="\n")
        extended_writer = csv.writer(
            extended_handle, delimiter="\t", lineterminator="\n"
        )
        final_writer.writerow(final_header)
        extended_writer.writerow(extended_header)

        for genotype in genotype_cols:
            population = popmap.get(genotype, "NA")
            roh_bed = roh_dir / "{}.bed".format(genotype)
            has_roh = roh_bed.exists() and roh_bed.stat().st_size > 0

            if has_roh:
                roh_ids = get_roh_variant_ids(variants_bed, roh_bed)
            else:
                roh_ids = set()

            if genotype in individual_lengths:
                roh_cds_bp, nonroh_cds_bp = individual_lengths[genotype]
            else:
                roh_cds_bp = 0
                nonroh_cds_bp = total_cds_bp

            roh_dosage_sum = 0.0
            nonroh_dosage_sum = 0.0
            roh_nonmissing = 0
            nonroh_nonmissing = 0
            roh_n_variants = len(roh_ids)
            nonroh_n_variants = n_variants - roh_n_variants

            for variant_id in all_variant_ids:
                dosage = variants[variant_id][genotype]
                if variant_id in roh_ids:
                    if dosage is not None:
                        roh_dosage_sum += dosage
                        roh_nonmissing += 1
                else:
                    if dosage is not None:
                        nonroh_dosage_sum += dosage
                        nonroh_nonmissing += 1

            if roh_cds_bp > 0:
                roh_rate = roh_dosage_sum / roh_cds_bp * 100000
                roh_rate_str = str(roh_rate)
            else:
                roh_rate = None
                roh_rate_str = "NA"

            if nonroh_cds_bp > 0:
                nonroh_rate = nonroh_dosage_sum / nonroh_cds_bp * 100000
                nonroh_rate_str = str(nonroh_rate)
            else:
                nonroh_rate = None
                nonroh_rate_str = "NA"

            final_writer.writerow(
                [genotype, population, roh_rate_str, nonroh_rate_str, roh_cds_bp, nonroh_cds_bp]
            )
            extended_writer.writerow(
                [
                    genotype,
                    population,
                    roh_dosage_sum,
                    nonroh_dosage_sum,
                    roh_cds_bp,
                    nonroh_cds_bp,
                    roh_n_variants,
                    nonroh_n_variants,
                    n_variants,
                    roh_nonmissing,
                    nonroh_nonmissing,
                    roh_rate_str,
                    nonroh_rate_str,
                ]
            )

    with final_output.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) != 6:
            fail(
                "Final output validation failed: file does not have exactly six columns."
            )

    print("CDS load variants: {}".format(n_variants))
    print("total_CDS_bp: {}".format(total_cds_bp))
    print("Wrote {}".format(final_output))
    print("Wrote {}".format(extended_qc))


if __name__ == "__main__":
    main()
