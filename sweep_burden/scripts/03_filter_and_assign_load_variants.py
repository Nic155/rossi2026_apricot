#!/usr/bin/env python3
"""
Filter PROVEAN/HIGH load variants and assign them to sweep or control CDS.

The large genotype matrix is streamed with Python's standard csv module.
Temporary BED intervals are used with bedtools intersect to classify variants.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path


WORKDIR = Path("~/work/sweep_gl").expanduser()
DEFAULT_INPUT = WORKDIR / "result3-provean_forGL+score+class.C1_W2.subset.tsv"
META_COLS = ["CHR", "POS", "AA", "IMPACT", "EFFECT", "PROVEAN_SCORE", "DEL_CLASS"]
LOAD_IMPACTS = {"PROVEAN", "HIGH"}
VARIANT_ID_COL = "__variant_id"


def fail(message):
    print("ERROR: {}".format(message), file=sys.stderr)
    sys.exit(1)


def require_file(path, label):
    if not path.exists() or path.stat().st_size == 0:
        fail("{} does not exist or is empty: {}".format(label, path))


def require_columns(columns):
    missing = [col for col in META_COLS if col not in columns]
    if missing:
        fail("Missing required metadata column(s): {}".format(", ".join(missing)))


def run_intersect(input_bed, region_bed, output_bed):
    with output_bed.open("w") as out_handle:
        subprocess.run(
            ["bedtools", "intersect", "-u", "-a", str(input_bed), "-b", str(region_bed)],
            check=True,
            stdout=out_handle,
        )


def read_variant_ids(path):
    ids = set()
    if not path.exists() or path.stat().st_size == 0:
        return ids
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                fail("Malformed BED line in {}: {}".format(path, line.rstrip()))
            ids.add(fields[3])
    return ids


def write_assigned_tables(filtered_table, original_header, sweep_ids, control_ids, sweep_output, control_output):
    sweep_rows = 0
    control_rows = 0

    with filtered_table.open(newline="") as in_handle, \
            sweep_output.open("w", newline="") as sweep_handle, \
            control_output.open("w", newline="") as control_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        if reader.fieldnames is None or VARIANT_ID_COL not in reader.fieldnames:
            fail("Temporary filtered table is missing {}.".format(VARIANT_ID_COL))

        sweep_writer = csv.DictWriter(
            sweep_handle,
            fieldnames=original_header,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        control_writer = csv.DictWriter(
            control_handle,
            fieldnames=original_header,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        sweep_writer.writeheader()
        control_writer.writeheader()

        for row in reader:
            variant_id = row[VARIANT_ID_COL]
            if variant_id in sweep_ids:
                sweep_writer.writerow(row)
                sweep_rows += 1
            elif variant_id in control_ids:
                control_writer.writerow(row)
                control_rows += 1

    return sweep_rows, control_rows


def write_qc(qc_path, metrics):
    with qc_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "count"])
        for key, value in metrics:
            writer.writerow([key, value])


def main():
    parser = argparse.ArgumentParser(
        description="Filter deleterious variants and assign them to sweep/control CDS."
    )
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Genetic-load TSV")
    parser.add_argument("--sweep-cds", default="results/regions/sweep_cds.bed")
    parser.add_argument("--control-cds", default="results/regions/control_cds.bed")
    parser.add_argument("--sweep-output", default="results/variants/load_variants_sweep.tsv")
    parser.add_argument("--control-output", default="results/variants/load_variants_control.tsv")
    parser.add_argument("--qc", default="results/qc/load_variant_counts.tsv")
    args = parser.parse_args()

    os.chdir(str(WORKDIR))

    input_path = Path(args.input).expanduser()
    sweep_cds = Path(args.sweep_cds)
    control_cds = Path(args.control_cds)
    sweep_output = Path(args.sweep_output)
    control_output = Path(args.control_output)
    qc_path = Path(args.qc)

    require_file(input_path, "Genetic-load TSV")
    require_file(sweep_cds, "Sweep CDS BED")
    require_file(control_cds, "Control CDS BED")
    if shutil.which("bedtools") is None:
        fail("bedtools was not found in PATH. Load bedtools before running this step.")

    Path("results/variants").mkdir(parents=True, exist_ok=True)
    Path("results/qc").mkdir(parents=True, exist_ok=True)
    tmp_dir = Path("results/variants/tmp_assign")
    tmp_dir.mkdir(parents=True, exist_ok=True)

    filtered_table = tmp_dir / "load_variants_filtered_with_id.tsv"
    variant_bed = tmp_dir / "load_variants_filtered.bed"
    sweep_hits_bed = tmp_dir / "load_variants_sweep_hits.bed"
    control_hits_bed = tmp_dir / "load_variants_control_hits.bed"

    for path in [filtered_table, variant_bed, sweep_hits_bed, control_hits_bed]:
        if path.exists():
            path.unlink()

    print("Reading and filtering load variants from {}".format(input_path))
    total_variants = 0
    load_variants = 0
    next_variant_id = 1

    with input_path.open(newline="") as in_handle, \
            filtered_table.open("w", newline="") as filtered_handle, \
            variant_bed.open("w", newline="") as bed_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        if reader.fieldnames is None:
            fail("Genetic-load TSV has no header.")
        header = list(reader.fieldnames)
        require_columns(header)
        genotype_cols = [col for col in header if col not in META_COLS]
        if not genotype_cols:
            fail("No genotype columns were found after metadata columns.")

        filtered_header = header + [VARIANT_ID_COL]
        filtered_writer = csv.DictWriter(
            filtered_handle,
            fieldnames=filtered_header,
            delimiter="\t",
            lineterminator="\n",
        )
        filtered_writer.writeheader()
        bed_writer = csv.writer(bed_handle, delimiter="\t", lineterminator="\n")

        for line_number, row in enumerate(reader, start=2):
            total_variants += 1
            if row.get("IMPACT") not in LOAD_IMPACTS:
                continue

            try:
                pos = int(float(row["POS"]))
            except ValueError:
                fail("A PROVEAN/HIGH variant has non-numeric POS at line {}".format(line_number))
            if pos < 1:
                fail("A PROVEAN/HIGH variant has POS < 1 at line {}".format(line_number))

            variant_id = "var_{}".format(next_variant_id)
            next_variant_id += 1
            load_variants += 1
            row[VARIANT_ID_COL] = variant_id

            filtered_writer.writerow(row)
            bed_writer.writerow([row["CHR"], pos - 1, pos, variant_id])

    if load_variants == 0:
        fail('No variants had IMPACT == "PROVEAN" or IMPACT == "HIGH".')

    print("Total variants in input: {}".format(total_variants))
    print("PROVEAN/HIGH variants: {}".format(load_variants))
    print("Assigning filtered variants to sweep/control CDS with bedtools intersect")
    run_intersect(variant_bed, sweep_cds, sweep_hits_bed)
    run_intersect(variant_bed, control_cds, control_hits_bed)

    sweep_ids = read_variant_ids(sweep_hits_bed)
    control_ids = read_variant_ids(control_hits_bed)
    both = sweep_ids.intersection(control_ids)
    if both:
        examples = ", ".join(sorted(list(both))[:10])
        fail(
            "Some variants overlap both sweep and control CDS, which should be "
            "disjoint. Example variant IDs: {}".format(examples)
        )

    any_cds_ids = sweep_ids.union(control_ids)
    discarded = load_variants - len(any_cds_ids)

    sweep_rows, control_rows = write_assigned_tables(
        filtered_table,
        header,
        sweep_ids,
        control_ids,
        sweep_output,
        control_output,
    )

    if sweep_rows != len(sweep_ids):
        fail("Mismatch between sweep BED hits and written sweep variant rows.")
    if control_rows != len(control_ids):
        fail("Mismatch between control BED hits and written control variant rows.")

    write_qc(
        qc_path,
        [
            ("total_variants_in_input", total_variants),
            ("PROVEAN_HIGH_variants", load_variants),
            ("PROVEAN_HIGH_variants_overlapping_any_CDS", len(any_cds_ids)),
            ("variants_assigned_to_sweep", len(sweep_ids)),
            ("variants_assigned_to_control", len(control_ids)),
            ("variants_discarded_because_outside_CDS", discarded),
        ],
    )

    print("PROVEAN/HIGH variants overlapping any CDS: {}".format(len(any_cds_ids)))
    print("Variants assigned to sweep: {}".format(len(sweep_ids)))
    print("Variants assigned to control: {}".format(len(control_ids)))
    print("Variants discarded because outside CDS: {}".format(discarded))
    print("Wrote {}".format(sweep_output))
    print("Wrote {}".format(control_output))
    print("Wrote {}".format(qc_path))


if __name__ == "__main__":
    main()
