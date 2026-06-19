#!/usr/bin/env python3
"""
Create merged candidate sweep regions from pi-ratio windows.

Selected windows are 1-based inclusive in the input CSV and are converted to
0-based half-open BED intervals before bedtools merge.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path


WORKDIR = Path("~/work/sweep_gl").expanduser()
DEFAULT_INPUT = WORKDIR / "pi_ratio_W2_C1.csv"
DEFAULT_THRESHOLD_COL = "95% threshold"
MERGE_DISTANCE_BP = 10000


def fail(message):
    print("ERROR: {}".format(message), file=sys.stderr)
    sys.exit(1)


def normalize_chr(value):
    chrom = str(value).strip()
    if not chrom or chrom.lower() == "nan":
        fail("Encountered an empty chromosome value in the pi-ratio file.")
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


def is_selected(value):
    text = str(value).strip()
    return text != "" and text.upper() not in {"NA", "NAN", "NULL"}


def run_bedtools_merge(selected_bed, merged_bed):
    if shutil.which("bedtools") is None:
        fail("bedtools was not found in PATH. Load bedtools before running this step.")

    with merged_bed.open("w") as out_handle:
        subprocess.run(
            [
                "bedtools",
                "merge",
                "-d",
                str(MERGE_DISTANCE_BP),
                "-i",
                str(selected_bed),
            ],
            check=True,
            stdout=out_handle,
        )


def bed_total_bp(path):
    total = 0
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            total += int(fields[2]) - int(fields[1])
    return total


def bed_count(path):
    with path.open() as handle:
        return sum(1 for line in handle if line.strip())


def main():
    parser = argparse.ArgumentParser(
        description="Select pi-ratio sweep windows and merge them with bedtools."
    )
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="pi_ratio_W2_C1.csv")
    parser.add_argument(
        "--threshold-col",
        default=DEFAULT_THRESHOLD_COL,
        help='Threshold-indicator column to use. Default: "95%% threshold"',
    )
    parser.add_argument(
        "--selected-bed",
        default="results/regions/sweep_windows_selected.bed",
        help="Output BED for selected 10 kb windows.",
    )
    parser.add_argument(
        "--merged-bed",
        default="results/regions/sweep_regions_merged.bed",
        help="Output BED for merged sweep regions.",
    )
    parser.add_argument(
        "--qc",
        default="results/qc/sweep_region_qc.tsv",
        help="Output QC table.",
    )
    args = parser.parse_args()

    os.chdir(str(WORKDIR))
    input_path = Path(args.input).expanduser()
    selected_bed = Path(args.selected_bed)
    merged_bed = Path(args.merged_bed)
    qc_path = Path(args.qc)

    if not input_path.exists():
        fail("Input file does not exist: {}".format(input_path))

    selected_bed.parent.mkdir(parents=True, exist_ok=True)
    merged_bed.parent.mkdir(parents=True, exist_ok=True)
    qc_path.parent.mkdir(parents=True, exist_ok=True)

    print("Reading pi-ratio windows: {}".format(input_path))
    intervals = []
    with input_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        if reader.fieldnames is None:
            fail("pi-ratio file has no header.")
        reader.fieldnames = [field.lstrip("\ufeff") for field in reader.fieldnames]
        required = ["Chr", "pos_start", "pos_end", args.threshold_col]
        missing = [col for col in required if col not in reader.fieldnames]
        if missing:
            fail("Missing required column(s) in pi-ratio file: {}".format(", ".join(missing)))

        data_started = False
        stopped_at_separator = False
        for line_number, row in enumerate(reader, start=2):
            row_values = [str(row.get(field, "")).strip() for field in reader.fieldnames]
            if all(value == "" for value in row_values):
                if data_started:
                    stopped_at_separator = True
                    print(
                        "Stopping at blank separator line {} after first pi-ratio table block.".format(
                            line_number
                        )
                    )
                    break
                continue

            data_started = True
            if not is_selected(row.get(args.threshold_col, "")):
                continue
            try:
                start = int(float(row["pos_start"])) - 1
                end = int(float(row["pos_end"]))
            except ValueError:
                fail("Non-numeric pos_start or pos_end at line {}".format(line_number))
            if start < 0 or end <= start:
                fail("Invalid BED coordinates generated at line {}".format(line_number))
            chrom = normalize_chr(row["Chr"])
            intervals.append((chrom, start, end))

    if not intervals:
        fail("No sweep windows passed threshold column: {}".format(args.threshold_col))

    if not stopped_at_separator:
        print("Read to end of pi-ratio file without finding a blank separator after the first table block.")

    intervals.sort(key=lambda item: (chrom_sort_key(item[0]), item[1], item[2]))
    with selected_bed.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(intervals)

    selected_count = len(intervals)
    selected_bp = sum(end - start for _chrom, start, end in intervals)
    print("Selected 10 kb windows: {}".format(selected_count))
    print("Total selected-window bp: {}".format(selected_bp))

    print("Merging selected windows with bedtools merge -d {}".format(MERGE_DISTANCE_BP))
    run_bedtools_merge(selected_bed, merged_bed)

    merged_count = bed_count(merged_bed)
    merged_bp = bed_total_bp(merged_bed)
    print("Merged sweep regions: {}".format(merged_count))
    print("Total merged sweep-region bp: {}".format(merged_bp))

    with qc_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["threshold_column", args.threshold_col])
        writer.writerow(["merge_distance_bp", MERGE_DISTANCE_BP])
        writer.writerow(["selected_10kb_windows", selected_count])
        writer.writerow(["total_selected_window_bp", selected_bp])
        writer.writerow(["merged_sweep_regions", merged_count])
        writer.writerow(["total_merged_sweep_region_bp", merged_bp])

    print("Wrote {}".format(selected_bed))
    print("Wrote {}".format(merged_bed))
    print("Wrote {}".format(qc_path))


if __name__ == "__main__":
    main()
