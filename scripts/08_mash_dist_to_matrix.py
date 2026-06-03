#!/usr/bin/env python3
"""Convert Mash pairwise distances to a square distance matrix."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from pipeline_common import parse_fasta, write_tsv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--mash-dist", required=True, type=Path)
    parser.add_argument("--out-matrix", required=True, type=Path)
    parser.add_argument("--out-long", required=True, type=Path)
    return parser.parse_args()


def load_tip_ids(fasta: Path) -> list[str]:
    tips = [header.split()[0] for header, _ in parse_fasta(fasta)]
    if len(tips) != len(set(tips)):
        raise ValueError(f"Duplicate FASTA IDs in {fasta}")
    return tips


def clean_mash_id(raw_id: str) -> str:
    return Path(raw_id).name.split()[0]


def main() -> None:
    args = parse_args()
    tips = load_tip_ids(args.fasta)
    tip_set = set(tips)
    distances = {(tip, tip): 0.0 for tip in tips}
    long_rows: list[dict] = []

    with args.mash_dist.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                raise ValueError(
                    f"Expected at least 5 Mash columns in {args.mash_dist}:{line_number}"
                )

            query = clean_mash_id(parts[0])
            reference = clean_mash_id(parts[1])
            if query not in tip_set or reference not in tip_set:
                raise ValueError(
                    f"Unexpected Mash ID at line {line_number}: {query!r}, {reference!r}"
                )

            distance = float(parts[2])
            if not math.isfinite(distance):
                raise ValueError(f"Non-finite Mash distance at line {line_number}: {parts[2]}")

            p_value = parts[3]
            shared_hashes = parts[4]
            distances[(query, reference)] = distance
            distances[(reference, query)] = distance
            long_rows.append(
                {
                    "query": query,
                    "reference": reference,
                    "distance": f"{distance:.12g}",
                    "p_value": p_value,
                    "shared_hashes": shared_hashes,
                }
            )

    missing_pairs = [
        f"{left}-{right}"
        for left in tips
        for right in tips
        if (left, right) not in distances
    ]
    if missing_pairs:
        preview = ", ".join(missing_pairs[:10])
        raise ValueError(
            f"Mash output did not include all pairwise distances; missing {len(missing_pairs)} "
            f"matrix cells. First missing: {preview}"
        )

    matrix_rows = []
    for tip in tips:
        row = {"tip": tip}
        for other in tips:
            row[other] = f"{distances[(tip, other)]:.12g}"
        matrix_rows.append(row)

    write_tsv(args.out_matrix, matrix_rows, ["tip", *tips])
    write_tsv(args.out_long, long_rows, ["query", "reference", "distance", "p_value", "shared_hashes"])
    print(f"Wrote {args.out_matrix}")
    print(f"Wrote {args.out_long}")


if __name__ == "__main__":
    main()
