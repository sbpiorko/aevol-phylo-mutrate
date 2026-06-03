#!/usr/bin/env python3
"""Summarize a multiple sequence alignment FASTA."""

from __future__ import annotations

import argparse
from pathlib import Path

from pipeline_common import parse_fasta, write_tsv


VALID_ALIGNMENT = set("ACGTN-")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args()


def is_parsimony_informative(column: list[str]) -> bool:
    counts: dict[str, int] = {}
    for char in column:
        if char in {"-", "N"}:
            continue
        counts[char] = counts.get(char, 0) + 1
    return sum(1 for count in counts.values() if count >= 2) >= 2


def main() -> None:
    args = parse_args()
    records = parse_fasta(args.alignment)
    if not records:
        raise ValueError(f"No alignment records found in {args.alignment}")

    lengths = [len(sequence) for _, sequence in records]
    if len(set(lengths)) != 1:
        raise ValueError(f"Alignment sequences do not all have the same length: {args.alignment}")

    invalid = sorted(
        {
            char
            for _, sequence in records
            for char in sequence
            if char not in VALID_ALIGNMENT
        }
    )
    if invalid:
        raise ValueError(f"Invalid alignment character(s): {''.join(invalid)}")

    alignment_length = lengths[0]
    variable_sites = 0
    parsimony_informative_sites = 0
    gap_columns = 0
    total_gaps = 0
    total_n = 0

    for i in range(alignment_length):
        column = [sequence[i] for _, sequence in records]
        states_without_gaps = {char for char in column if char not in {"-", "N"}}
        if len(states_without_gaps) > 1:
            variable_sites += 1
        if is_parsimony_informative(column):
            parsimony_informative_sites += 1
        if "-" in column:
            gap_columns += 1
        total_gaps += column.count("-")
        total_n += column.count("N")

    total_cells = alignment_length * len(records)
    row = {
        "alignment": args.alignment.as_posix(),
        "records": len(records),
        "alignment_length": alignment_length,
        "variable_sites": variable_sites,
        "parsimony_informative_sites": parsimony_informative_sites,
        "gap_columns": gap_columns,
        "gap_fraction": f"{(total_gaps / total_cells):.6f}" if total_cells else "0.000000",
        "ambiguous_n_fraction": f"{(total_n / total_cells):.6f}" if total_cells else "0.000000",
    }
    write_tsv(
        args.out,
        [row],
        [
            "alignment",
            "records",
            "alignment_length",
            "variable_sites",
            "parsimony_informative_sites",
            "gap_columns",
            "gap_fraction",
            "ambiguous_n_fraction",
        ],
    )
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
