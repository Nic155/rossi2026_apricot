#!/usr/bin/env python3
"""
Filter PROVEAN/HIGH load variants and restrict to CDS regions.

Variants with IMPACT == "PROVEAN" or IMPACT == "HIGH" are selected.
The selected variants are intersected with the reduced CDS BED using bedtools
to remove variants that do not overlap any CDS interval.

Outputs:
  results/variants/load_variants_cds.tsv -- filtered CDS variants (original columns)
  results/variants/load_variants_cds.bed -- BED with variant_id in col 4 (used by step 04)
  results/qc/load_variant_counts.tsv     -- QC counts
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path


WORKDIR = Path("~/work/ROH/roh_gl").expanduser()
DEFAULT_INPUT = Path(
    "~/work/ROH/result3-provean_forGL+score+class.C1_W2.subset.tsv"
).expanduser()
META_COLS = ["CHR", "POS", "AA", "IMPACT", "EFFECT", "PROVEAN_SCORE", "DEL_CLASS"]
LOAD_IMPACTS = {"PROVEAN", "HIGH"}
VARIANT_ID_COL = "__variant_id"


def fail(message):
    print("ERROR: {}".format(message), file=sys.stderr)
    sys.exit(1)


def require_file(path, label):
    if not path.exists() or path.stat().st_size == 0:
        fail("{} does not exist or is empty: {}".format(label, path))


def run_intersect(input_bed, region_bed, output_bed):
    with output_bed.open("w") as out_handle:
        subprocess.run(
            ["bedtools", "intersect", "-u", "-a", str(input_bed), "-b", str(region_bed)],
            check=True,
            stdout=out_handle,
        )


def read_hit_ids(path):
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


def write_qc(qc_path, metrics):
    with qc_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "count"])
        for key, value in metrics:
            writer.writerow([key, value])


def main():
    parser = argparse.ArgumentParser(
        description="Filter deleterious variants and restrict to CDS regions."
    )
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Genetic-load TSV")
    parser.add_argument("--cds", default="results/regions/cds_reduced.bed")
    parser.add_argument("--output", default="results/variants/load_variants_cds.tsv")
    parser.add_argument("--output-bed", default="results/variants/load_variants_cds.bed")
    parser.add_argument("--qc", default="results/qc/load_variant_counts.tsv")
    args = parser.parse_args()

    os.chdir(str(WORKDIR))

    input_path = Path(args.input).expanduser()
    cds_path = Path(args.cds)
    output_path = Path(args.output)
    output_bed = Path(args.output_bed)
    qc_path = Path(args.qc)

    require_file(input_path, "Genetic-load TSV")
    require_file(cds_path, "Reduced CDS BED")
    if shutil.which("bedtools") is None:
        fail("bedtools was not found in PATH. Load bedtools before running this step.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir = Path("results/variants/tmp_filter")
    tmp_dir.mkdir(parents=True, exist_ok=True)

    tmp_all_bed = tmp_dir / "load_variants_all.bed"
    tmp_all_tsv = tmp_dir / "load_variants_all_with_id.tsv"
    tmp_cds_hits = tmp_dir / "load_variants_cds_hits.bed"

    for path in [tmp_all_bed, tmp_all_tsv, tmp_cds_hits]:
        if path.exists():
            path.unlink()

    print("Reading and filtering load variants from {}".format(input_path))
    total_variants = 0
    load_variants = 0
    next_variant_id = 1

    with input_path.open(newline="") as in_handle, \
            tmp_all_tsv.open("w", newline="") as tsv_handle, \
            tmp_all_bed.open("w", newline="") as bed_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        if reader.fieldnames is None:
            fail("Genetic-load TSV has no header.")
        header = list(reader.fieldnames)
        missing = [col for col in META_COLS if col not in header]
        if missing:
            fail("Missing required metadata column(s): {}".format(", ".join(missing)))
        genotype_cols = [col for col in header if col not in META_COLS]
        if not genotype_cols:
            fail("No genotype columns were found after metadata columns.")

        tsv_writer = csv.DictWriter(
            tsv_handle,
            fieldnames=header + [VARIANT_ID_COL],
            delimiter="\t",
            lineterminator="\n",
        )
        tsv_writer.writeheader()
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

            tsv_writer.writerow(row)
            bed_writer.writerow([row["CHR"], pos - 1, pos, variant_id])

    if load_variants == 0:
        fail('No variants had IMPACT == "PROVEAN" or IMPACT == "HIGH".')

    print("Total variants in input: {}".format(total_variants))
    print("PROVEAN/HIGH variants: {}".format(load_variants))
    print("Intersecting filtered variants with CDS using bedtools intersect")
    run_intersect(tmp_all_bed, cds_path, tmp_cds_hits)

    cds_ids = read_hit_ids(tmp_cds_hits)
    discarded = load_variants - len(cds_ids)

    cds_rows = 0
    with tmp_all_tsv.open(newline="") as in_handle, \
            output_path.open("w", newline="") as out_tsv_handle, \
            output_bed.open("w", newline="") as out_bed_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        if reader.fieldnames is None or VARIANT_ID_COL not in reader.fieldnames:
            fail("Temporary TSV is missing {}.".format(VARIANT_ID_COL))
        tsv_writer = csv.DictWriter(
            out_tsv_handle,
            fieldnames=header,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        bed_writer = csv.writer(out_bed_handle, delimiter="\t", lineterminator="\n")
        tsv_writer.writeheader()

        for row in reader:
            variant_id = row[VARIANT_ID_COL]
            if variant_id not in cds_ids:
                continue
            tsv_writer.writerow(row)
            pos = int(float(row["POS"]))
            bed_writer.writerow([row["CHR"], pos - 1, pos, variant_id])
            cds_rows += 1

    if cds_rows != len(cds_ids):
        fail(
            "Mismatch between CDS BED hits ({}) and written TSV rows ({}).".format(
                len(cds_ids), cds_rows
            )
        )

    write_qc(
        qc_path,
        [
            ("total_variants_in_input", total_variants),
            ("PROVEAN_HIGH_variants", load_variants),
            ("PROVEAN_HIGH_variants_in_CDS", len(cds_ids)),
            ("variants_discarded_outside_CDS", discarded),
        ],
    )

    print("PROVEAN/HIGH variants in CDS: {}".format(len(cds_ids)))
    print("Variants discarded because outside CDS: {}".format(discarded))
    print("Wrote {}".format(output_path))
    print("Wrote {}".format(output_bed))
    print("Wrote {}".format(qc_path))


if __name__ == "__main__":
    main()
