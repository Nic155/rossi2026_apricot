#!/usr/bin/env python3
"""
Calculate per-genotype deleterious dosage per 100 kb CDS for sweep and control.
"""

import argparse
import csv
import os
import sys
from pathlib import Path


WORKDIR = Path("~/work/sweep_gl").expanduser()
META_COLS = ["CHR", "POS", "AA", "IMPACT", "EFFECT", "PROVEAN_SCORE", "DEL_CLASS"]
MISSING_VALUES = {"", "NA", "NAN", "NULL", "."}


def fail(message):
    print("ERROR: {}".format(message), file=sys.stderr)
    sys.exit(1)


def require_file(path, label):
    if not path.exists() or path.stat().st_size == 0:
        fail("{} does not exist or is empty: {}".format(label, path))


def read_cds_lengths(path):
    length_map = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail("CDS lengths file has no header.")
        required = {"region_class", "cds_length_bp"}
        missing = required.difference(set(reader.fieldnames))
        if missing:
            fail("CDS lengths file is missing column(s): {}".format(", ".join(sorted(missing))))
        for row in reader:
            length_map[row["region_class"]] = int(float(row["cds_length_bp"]))

    for key in ["sweep", "control"]:
        if key not in length_map:
            fail("CDS lengths file is missing region_class: {}".format(key))

    sweep_bp = length_map["sweep"]
    control_bp = length_map["control"]
    if sweep_bp <= 0:
        fail("sweep_CDS_bp is zero.")
    if control_bp <= 0:
        fail("control_CDS_bp is zero.")
    return sweep_bp, control_bp


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


def summarize_variant_table(path, genotype_cols):
    dosage_sum = {genotype: 0.0 for genotype in genotype_cols}
    nonmissing_calls = {genotype: 0 for genotype in genotype_cols}
    n_variants = 0

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail("{} has no header.".format(path))
        missing = [col for col in genotype_cols if col not in reader.fieldnames]
        if missing:
            fail("{} is missing genotype column(s): {}".format(path, ", ".join(missing)))

        for line_number, row in enumerate(reader, start=2):
            n_variants += 1
            for genotype in genotype_cols:
                dosage = parse_dosage(row.get(genotype, ""), path, line_number, genotype)
                if dosage is None:
                    continue
                dosage_sum[genotype] += dosage
                nonmissing_calls[genotype] += 1

    return n_variants, dosage_sum, nonmissing_calls


def main():
    parser = argparse.ArgumentParser(
        description="Calculate soybean-style per-genotype deleterious burden."
    )
    parser.add_argument("--sweep-variants", default="results/variants/load_variants_sweep.tsv")
    parser.add_argument("--control-variants", default="results/variants/load_variants_control.tsv")
    parser.add_argument("--cds-lengths", default="results/qc/cds_lengths.tsv")
    parser.add_argument("--final-output", default="results/final_genotype_sweep_control_burden.tsv")
    parser.add_argument("--extended-qc", default="results/qc/genotype_burden_extended.tsv")
    args = parser.parse_args()

    os.chdir(str(WORKDIR))

    sweep_variants = Path(args.sweep_variants)
    control_variants = Path(args.control_variants)
    cds_lengths = Path(args.cds_lengths)
    final_output = Path(args.final_output)
    extended_qc = Path(args.extended_qc)

    require_file(sweep_variants, "Sweep variant table")
    require_file(control_variants, "Control variant table")
    require_file(cds_lengths, "CDS lengths table")

    final_output.parent.mkdir(parents=True, exist_ok=True)
    extended_qc.parent.mkdir(parents=True, exist_ok=True)

    sweep_bp, control_bp = read_cds_lengths(cds_lengths)
    sweep_genotypes = get_genotype_columns(sweep_variants)
    control_genotypes = get_genotype_columns(control_variants)

    if sweep_genotypes != control_genotypes:
        fail("Sweep and control variant tables do not have identical genotype columns.")

    genotype_cols = sweep_genotypes
    print("Genotype count: {}".format(len(genotype_cols)))

    sweep_n_variants, sweep_sum, sweep_calls = summarize_variant_table(
        sweep_variants, genotype_cols
    )
    control_n_variants, control_sum, control_calls = summarize_variant_table(
        control_variants, genotype_cols
    )

    if sweep_n_variants == 0:
        fail("No sweep variants are present.")
    if control_n_variants == 0:
        fail("No control variants are present.")

    final_header = [
        "genotype",
        "sweep_deleterious_dosage_per_100kb_CDS",
        "control_deleterious_dosage_per_100kb_CDS",
    ]
    extended_header = [
        "genotype",
        "sweep_raw_dosage_sum",
        "control_raw_dosage_sum",
        "sweep_CDS_bp",
        "control_CDS_bp",
        "sweep_n_variants",
        "control_n_variants",
        "sweep_n_nonmissing_calls",
        "control_n_nonmissing_calls",
        "sweep_deleterious_dosage_per_100kb_CDS",
        "control_deleterious_dosage_per_100kb_CDS",
    ]

    with final_output.open("w", newline="") as final_handle, \
            extended_qc.open("w", newline="") as extended_handle:
        final_writer = csv.writer(final_handle, delimiter="\t", lineterminator="\n")
        extended_writer = csv.writer(extended_handle, delimiter="\t", lineterminator="\n")
        final_writer.writerow(final_header)
        extended_writer.writerow(extended_header)

        for genotype in genotype_cols:
            sweep_rate = sweep_sum[genotype] / sweep_bp * 100000
            control_rate = control_sum[genotype] / control_bp * 100000
            final_writer.writerow([genotype, sweep_rate, control_rate])
            extended_writer.writerow(
                [
                    genotype,
                    sweep_sum[genotype],
                    control_sum[genotype],
                    sweep_bp,
                    control_bp,
                    sweep_n_variants,
                    control_n_variants,
                    sweep_calls[genotype],
                    control_calls[genotype],
                    sweep_rate,
                    control_rate,
                ]
            )

    with final_output.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) != 3:
            fail("Final output validation failed: file does not have exactly three columns.")

    print("Sweep variants: {}".format(sweep_n_variants))
    print("Control variants: {}".format(control_n_variants))
    print("sweep_CDS_bp: {}".format(sweep_bp))
    print("control_CDS_bp: {}".format(control_bp))
    print("Wrote {}".format(final_output))
    print("Wrote {}".format(extended_qc))


if __name__ == "__main__":
    main()
