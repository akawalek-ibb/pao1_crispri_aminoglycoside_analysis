#!/usr/bin/env python3
"""
Combine per-gene CIAO predictions, preserving original columns and raw scores.
Prerequisite: 02_run_CIAO_batch.py or existing per-gene prediction folders.
Input: <root>/<locus_tag>/gRNAs.csv (usually tab-separated).
Output: one combined CSV. Requires pandas; no prediction model is loaded.
Run from the project root:
  python guide_efficacy_prediction/ciao/03_concatenate_CIAO_gRNA_predictions.py \
    --root <results> --id-column locus_tag --strict-files --output <combined.csv>
For PAO1 benchmarks, retain the default PAO1_ID column and optional --strict-pao1.
Temporary .*.running directories are excluded. Duplicate guide rows are reported
and retained. Use only after checking that the batch prediction run completed.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Optional

import pandas as pd


PAO1_PATTERN = re.compile(r"(?i)(PA\d{4})")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Recursively find gRNAs.csv files, add an identifier from each "
            "parent folder, and concatenate all predictions into one CSV."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="Parent folder to search recursively. Default: current folder.",
    )
    parser.add_argument(
        "--filename",
        default="gRNAs.csv",
        help="Prediction filename to search for. Default: gRNAs.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("all_PAO1_gRNA_predictions.csv"),
        help=(
            "Combined output CSV. Relative paths are resolved inside --root. "
            "Default: all_PAO1_gRNA_predictions.csv"
        ),
    )
    parser.add_argument(
        "--id-column", default="PAO1_ID",
        help="Output identifier column: keep PAO1_ID for PAO1 benchmarks, or use locus_tag for other strains.",
    )
    parser.add_argument(
        "--strict-pao1",
        action="store_true",
        help=(
            "Stop when a parent-folder name does not contain a PAO1 locus tag "
            "matching PA followed by four digits."
        ),
    )
    parser.add_argument(
        "--strict-files",
        action="store_true",
        help="Stop at the first unreadable or malformed prediction file.",
    )
    return parser.parse_args()


def detect_delimiter(path: Path) -> str:
    """Detect tab/comma/semicolon delimiter from the beginning of a file."""
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        sample = handle.read(8192)

    if not sample.strip():
        raise ValueError("file is empty")

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
        return dialect.delimiter
    except csv.Error:
        # CIAO output is normally tab-delimited even though its extension is .csv.
        first_line = sample.splitlines()[0]
        return "\t" if first_line.count("\t") >= first_line.count(",") else ","


def extract_pao1_id(folder_name: str, strict: bool) -> str:
    """Extract PA#### from a folder name, or retain the folder name as fallback."""
    match = PAO1_PATTERN.search(folder_name)
    if match:
        return match.group(1).upper()

    if strict:
        raise ValueError(
            f"parent folder '{folder_name}' does not contain a PAO1_ID such as PA0001"
        )

    # Fallback preserves data from unexpectedly named folders while making the
    # issue visible in the console warning.
    return folder_name


def read_prediction_file(path: Path, strict_pao1: bool, id_column: str = "PAO1_ID") -> pd.DataFrame:
    delimiter = detect_delimiter(path)
    table = pd.read_csv(
        path,
        sep=delimiter,
        dtype=str,
        encoding="utf-8-sig",
        keep_default_na=False,
    )

    if table.empty:
        raise ValueError("file contains a header but no guide rows")

    table.columns = [str(column).strip() for column in table.columns]

    required = {"gRNA_ID", "SequenceID"}
    missing = sorted(required.difference(table.columns))
    if missing:
        raise ValueError("missing required column(s): " + ", ".join(missing))

    parent_folder = path.parent.name
    pao1_id = extract_pao1_id(parent_folder, strict=strict_pao1) if id_column == "PAO1_ID" else parent_folder

    # Add identifiers at the start while retaining the original SequenceID.
    table.insert(0, "gene_name", table["SequenceID"].astype(str))
    table.insert(0, id_column, pao1_id)

    return table


def natural_pao1_key(value: str) -> tuple[int, str]:
    match = re.fullmatch(r"PA(\d{4})", str(value).upper())
    if match:
        return int(match.group(1)), str(value)
    return 10**9, str(value)


def main() -> int:
    args = parse_args()
    if args.id_column not in ("PAO1_ID", "locus_tag"):
        print("ERROR: --id-column must be PAO1_ID or locus_tag", file=sys.stderr)
        return 2
    if args.strict_pao1 and args.id_column != "PAO1_ID":
        print("ERROR: --strict-pao1 requires --id-column PAO1_ID", file=sys.stderr)
        return 2
    root = args.root.expanduser().resolve()

    if not root.is_dir():
        print(f"ERROR: root folder does not exist or is not a directory: {root}", file=sys.stderr)
        return 2

    output = args.output.expanduser()
    if not output.is_absolute():
        output = root / output
    output = output.resolve()

    prediction_files = sorted(
        (
            path
            for path in root.rglob(args.filename)
            if path.is_file() and path.resolve() != output
            and not any(part.startswith(".") and part.endswith(".running")
                        for part in path.relative_to(root).parts[:-1])
        ),
        key=lambda path: str(path.relative_to(root)),
    )

    if not prediction_files:
        print(
            f"ERROR: no files named '{args.filename}' were found below {root}",
            file=sys.stderr,
        )
        return 1

    tables: list[pd.DataFrame] = []
    failures: list[tuple[Path, str]] = []
    fallback_folders: list[str] = []

    for index, path in enumerate(prediction_files, start=1):
        relative_path = path.relative_to(root)
        try:
            table = read_prediction_file(path, strict_pao1=args.strict_pao1, id_column=args.id_column)
            if args.id_column == "PAO1_ID" and not PAO1_PATTERN.search(path.parent.name):
                fallback_folders.append(path.parent.name)
            tables.append(table)
            print(
                f"[{index}/{len(prediction_files)}] {relative_path}: "
                f"{len(table):,} guides"
            )
        except Exception as exc:  # report the exact file and continue unless strict
            failures.append((relative_path, str(exc)))
            print(f"WARNING: skipped {relative_path}: {exc}", file=sys.stderr)
            if args.strict_files or args.strict_pao1:
                return 1

    if not tables:
        print("ERROR: none of the discovered files could be read.", file=sys.stderr)
        return 1

    combined = pd.concat(tables, ignore_index=True, sort=False)

    # Deterministic ordering: locus tag, gene name, then activity score
    # where the latter is available. Numeric conversion is temporary only;
    # original text values remain unchanged in the output.
    combined["__pao1_order"] = combined[args.id_column].map(
        lambda value: natural_pao1_key(value)[0]
    )
    sort_columns = ["__pao1_order", args.id_column, "gene_name"]

    if "Activity score" in combined.columns:
        combined["__activity_order"] = pd.to_numeric(
            combined["Activity score"], errors="coerce"
        )
        sort_columns.append("__activity_order")

    combined = combined.sort_values(sort_columns, kind="stable", na_position="last")
    combined = combined.drop(
        columns=[column for column in ("__pao1_order", "__activity_order") if column in combined]
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output, index=False, encoding="utf-8")

    unique_genes = combined[args.id_column].nunique(dropna=True)
    duplicate_guides = 0
    if "gRNA sequence (20nt)" in combined.columns:
        duplicate_guides = int(
            combined.duplicated(
                subset=[args.id_column, "gRNA sequence (20nt)"],
                keep=False,
            ).sum()
        )

    print("\nCompleted")
    print(f"Files discovered:       {len(prediction_files):,}")
    print(f"Files concatenated:     {len(tables):,}")
    print(f"Files skipped:          {len(failures):,}")
    print(f"Unique {args.id_column} values:  {unique_genes:,}")
    print(f"Total guide rows:       {len(combined):,}")
    print(f"Duplicate guide rows*:  {duplicate_guides:,}")
    print(f"Output:                  {output}")
    print(f"*Same {args.id_column} and 20-nt guide sequence; rows are retained.")

    if fallback_folders:
        unique_fallbacks = sorted(set(fallback_folders))
        preview = ", ".join(unique_fallbacks[:10])
        suffix = " ..." if len(unique_fallbacks) > 10 else ""
        print(
            "WARNING: PA#### was not found in these folder names; the complete "
            f"folder name was used as PAO1_ID: {preview}{suffix}",
            file=sys.stderr,
        )

    if failures:
        print("\nSkipped files:", file=sys.stderr)
        for path, reason in failures:
            print(f"  - {path}: {reason}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
