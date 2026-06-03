#!/usr/bin/env python3
"""Standardize and validate FASTA inputs for phylogeny methods."""

from __future__ import annotations

import argparse
import hashlib
import re
from collections import Counter
from pathlib import Path

from pipeline_common import ensure_dir, parse_fasta, project_root, read_mutation_grid, write_fasta, write_tsv


VALID_DNA = set("ACGTN")
TIP_PATTERN = re.compile(r"^(S\d{2})(?:\b|[_\s.-].*)")


def parse_args() -> argparse.Namespace:
    root = project_root()
    default_conditions = list(read_mutation_grid(root).keys())
    parser = argparse.ArgumentParser(
        description="Prepare clean, validated FASTA files for phylogeny analysis."
    )
    parser.add_argument(
        "--conditions",
        nargs="+",
        default=default_conditions,
        help="Mutation-rate conditions to prepare. Default: all mutation-grid conditions.",
    )
    parser.add_argument(
        "--replicates",
        nargs="+",
        default=["rep_001"],
        help="Replicate labels to prepare. Default: rep_001.",
    )
    parser.add_argument("--expected-tips", type=int, default=16)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def clean_tip_id(header: str) -> str:
    """Return the Sxx tip label from an Aevol FASTA header."""
    header = header.strip()
    match = TIP_PATTERN.match(header)
    if match:
        return match.group(1)
    first_token = header.split()[0] if header.split() else header
    if re.fullmatch(r"S\d{2}", first_token):
        return first_token
    raise ValueError(f"Could not derive clean Sxx tip id from FASTA header: {header!r}")


def gc_fraction(sequence: str) -> float:
    if not sequence:
        return 0.0
    gc = sequence.count("G") + sequence.count("C")
    return gc / len(sequence)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def expected_tip_ids(expected_tips: int) -> list[str]:
    return [f"S{i:02d}" for i in range(1, expected_tips + 1)]


def validate_and_clean_records(
    records: list[tuple[str, str]],
    *,
    expected_tips: int,
    source: Path,
) -> tuple[list[tuple[str, str]], list[dict], dict]:
    clean_records: list[tuple[str, str]] = []
    sample_rows: list[dict] = []
    invalid_messages: list[str] = []
    seen_ids: list[str] = []
    sequence_hashes: list[str] = []

    for original_header, raw_sequence in records:
        tip_id = clean_tip_id(original_header)
        sequence = "".join(raw_sequence.split()).upper()
        seen_ids.append(tip_id)

        invalid_chars = sorted(set(sequence) - VALID_DNA)
        if invalid_chars:
            invalid_messages.append(
                f"{tip_id}: invalid DNA character(s): {''.join(invalid_chars)}"
            )
        if not sequence:
            invalid_messages.append(f"{tip_id}: empty sequence")

        seq_hash = sha256_text(sequence)
        sequence_hashes.append(seq_hash)
        clean_records.append((tip_id, sequence))
        sample_rows.append(
            {
                "tip": tip_id,
                "original_header": original_header,
                "length": len(sequence),
                "gc_fraction": f"{gc_fraction(sequence):.6f}",
                "ambiguous_n": sequence.count("N"),
                "sha256": seq_hash,
            }
        )

    expected_ids = expected_tip_ids(expected_tips)
    id_counts = Counter(seen_ids)
    duplicate_ids = sorted(tip_id for tip_id, count in id_counts.items() if count > 1)
    missing_ids = sorted(set(expected_ids) - set(seen_ids))
    extra_ids = sorted(set(seen_ids) - set(expected_ids))
    duplicate_sequence_count = sum(count - 1 for count in Counter(sequence_hashes).values() if count > 1)

    if len(records) != expected_tips:
        invalid_messages.append(
            f"expected {expected_tips} FASTA records but found {len(records)} in {source}"
        )
    if duplicate_ids:
        invalid_messages.append("duplicate tip id(s): " + ", ".join(duplicate_ids))
    if missing_ids:
        invalid_messages.append("missing expected tip id(s): " + ", ".join(missing_ids))
    if extra_ids:
        invalid_messages.append("unexpected tip id(s): " + ", ".join(extra_ids))

    clean_records = sorted(clean_records, key=lambda item: item[0])
    sample_rows = sorted(sample_rows, key=lambda row: row["tip"])

    lengths = [int(row["length"]) for row in sample_rows]
    summary = {
        "records": len(records),
        "min_length": min(lengths) if lengths else 0,
        "max_length": max(lengths) if lengths else 0,
        "mean_length": f"{(sum(lengths) / len(lengths)):.2f}" if lengths else "0.00",
        "duplicate_ids": len(duplicate_ids),
        "duplicate_sequences": duplicate_sequence_count,
        "invalid_records": len(invalid_messages),
        "status": "ok" if not invalid_messages else "error",
        "messages": "; ".join(invalid_messages),
    }
    return clean_records, sample_rows, summary


def prepare_one(
    *,
    root: Path,
    condition: str,
    replicate: str,
    expected_tips: int,
    force: bool,
) -> dict:
    raw_fasta = root / "data" / "raw_final_fastas" / condition / replicate / "all_tips.fa"
    clean_dir = ensure_dir(root / "data" / "clean_fastas" / condition / replicate)
    clean_fasta = clean_dir / "all_tips.clean.fa"
    manifest = clean_dir / "sample_manifest.tsv"

    if not raw_fasta.exists():
        raise FileNotFoundError(f"Missing raw FASTA: {raw_fasta}")
    if clean_fasta.exists() and manifest.exists() and not force:
        records = parse_fasta(clean_fasta)
        _, _, summary = validate_and_clean_records(
            records,
            expected_tips=expected_tips,
            source=clean_fasta,
        )
        return {
            "condition": condition,
            "replicate": replicate,
            **summary,
            "raw_fasta": raw_fasta.relative_to(root).as_posix(),
            "clean_fasta": clean_fasta.relative_to(root).as_posix(),
            "sample_manifest": manifest.relative_to(root).as_posix(),
            "action": "skipped_existing",
        }

    records = parse_fasta(raw_fasta)
    clean_records, sample_rows, summary = validate_and_clean_records(
        records,
        expected_tips=expected_tips,
        source=raw_fasta,
    )
    if summary["status"] != "ok":
        raise ValueError(f"{raw_fasta} failed validation: {summary['messages']}")

    write_fasta(clean_fasta, clean_records)
    write_tsv(
        manifest,
        sample_rows,
        ["tip", "original_header", "length", "gc_fraction", "ambiguous_n", "sha256"],
    )

    return {
        "condition": condition,
        "replicate": replicate,
        **summary,
        "raw_fasta": raw_fasta.relative_to(root).as_posix(),
        "clean_fasta": clean_fasta.relative_to(root).as_posix(),
        "sample_manifest": manifest.relative_to(root).as_posix(),
        "action": "written",
    }


def main() -> None:
    args = parse_args()
    root = project_root()
    summary_rows: list[dict] = []

    for condition in args.conditions:
        for replicate in args.replicates:
            row = prepare_one(
                root=root,
                condition=condition,
                replicate=replicate,
                expected_tips=args.expected_tips,
                force=args.force,
            )
            summary_rows.append(row)
            print(
                f"{condition} {replicate}: {row['status']} "
                f"({row['records']} records, mean length {row['mean_length']})"
            )

    summary_path = root / "results" / "simulation_qc" / "phylo_input_qc.tsv"
    write_tsv(
        summary_path,
        summary_rows,
        [
            "condition",
            "replicate",
            "records",
            "min_length",
            "max_length",
            "mean_length",
            "duplicate_ids",
            "duplicate_sequences",
            "invalid_records",
            "status",
            "messages",
            "raw_fasta",
            "clean_fasta",
            "sample_manifest",
            "action",
        ],
    )
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
