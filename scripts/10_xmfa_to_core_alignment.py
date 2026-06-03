#!/usr/bin/env python3
"""Convert progressiveMauve XMFA output to a concatenated core FASTA alignment."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from pipeline_common import write_fasta, write_tsv


SEQ_FILE_RE = re.compile(r"^#Sequence(\d+)File\s+(.+)$")
XMFA_HEADER_RE = re.compile(r"^>\s*(\d+):[^\s]*\s+([+-])\s*(.*)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xmfa", required=True, type=Path)
    parser.add_argument("--out-fasta", required=True, type=Path)
    parser.add_argument("--out-blocks", required=True, type=Path)
    parser.add_argument("--expected-tips", nargs="+", required=True)
    return parser.parse_args()


def tip_from_path(path_text: str) -> str:
    return Path(path_text.strip()).stem


def read_xmfa(path: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    sequence_files: dict[str, str] = {}
    blocks: list[dict[str, str]] = []
    current_block: dict[str, str] = {}
    current_seq: str | None = None
    chunks: list[str] = []

    def finish_record() -> None:
        nonlocal current_seq, chunks
        if current_seq is not None:
            if current_seq in current_block:
                raise ValueError(f"Duplicate sequence {current_seq} in one XMFA block")
            current_block[current_seq] = "".join(chunks).upper()
        current_seq = None
        chunks = []

    def finish_block() -> None:
        nonlocal current_block
        finish_record()
        if current_block:
            blocks.append(current_block)
        current_block = {}

    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        seq_file_match = SEQ_FILE_RE.match(line)
        if seq_file_match:
            sequence_files[seq_file_match.group(1)] = tip_from_path(seq_file_match.group(2))
            continue
        if line.startswith("#"):
            continue
        if line == "=":
            finish_block()
            continue
        if line.startswith(">"):
            finish_record()
            header_match = XMFA_HEADER_RE.match(line)
            if not header_match:
                raise ValueError(f"Could not parse XMFA sequence header: {line}")
            sequence_number = header_match.group(1)
            tip = sequence_files.get(sequence_number)
            if tip is None:
                trailing_path = header_match.group(3).strip()
                tip = tip_from_path(trailing_path) if trailing_path else f"Sequence{sequence_number}"
            current_seq = tip
            chunks = []
            continue
        if current_seq is None:
            raise ValueError(f"Sequence data encountered before XMFA header: {line[:40]}")
        chunks.append(line)

    finish_block()
    return sequence_files, blocks


def strip_all_gap_columns(block: dict[str, str], expected: list[str]) -> dict[str, str]:
    lengths = {len(block[tip]) for tip in expected}
    if len(lengths) != 1:
        raise ValueError("XMFA block has unequal aligned sequence lengths")
    block_length = lengths.pop()
    kept_columns = []
    for i in range(block_length):
        column = [block[tip][i] for tip in expected]
        if any(char != "-" for char in column):
            kept_columns.append(i)
    return {tip: "".join(block[tip][i] for i in kept_columns) for tip in expected}


def main() -> None:
    args = parse_args()
    expected = args.expected_tips
    _, blocks = read_xmfa(args.xmfa)
    missing_overall = []
    alignment_parts = {tip: [] for tip in expected}
    block_rows: list[dict] = []

    for block_index, block in enumerate(blocks, start=1):
        present = set(block)
        missing = [tip for tip in expected if tip not in present]
        extra = sorted(present - set(expected))
        used = False
        ungapped_length = 0
        aligned_length = 0
        if not missing and not extra:
            stripped = strip_all_gap_columns(block, expected)
            aligned_length = len(next(iter(stripped.values()))) if stripped else 0
            ungapped_length = sum(1 for i in range(aligned_length) if any(stripped[tip][i] != "-" for tip in expected))
            if aligned_length > 0:
                for tip in expected:
                    alignment_parts[tip].append(stripped[tip])
                used = True
        else:
            missing_overall.extend(missing)

        block_rows.append(
            {
                "block": block_index,
                "used_in_core": str(used).lower(),
                "present_tips": len(present.intersection(expected)),
                "missing_tips": ",".join(missing),
                "extra_tips": ",".join(extra),
                "aligned_length": aligned_length,
                "ungapped_nonempty_columns": ungapped_length,
            }
        )

    records = [(tip, "".join(alignment_parts[tip])) for tip in expected]
    lengths = {len(sequence) for _, sequence in records}
    if not records or lengths == {0}:
        raise ValueError(f"No core alignment columns were recovered from {args.xmfa}")
    if len(lengths) != 1:
        raise ValueError("Concatenated core alignment has unequal sequence lengths")
    all_gap_tips = [tip for tip, sequence in records if set(sequence) <= {"-"}]
    if all_gap_tips:
        raise ValueError(
            "Core alignment leaves these tips with only gaps: "
            + ", ".join(all_gap_tips)
        )

    write_fasta(args.out_fasta, records)
    write_tsv(
        args.out_blocks,
        block_rows,
        [
            "block",
            "used_in_core",
            "present_tips",
            "missing_tips",
            "extra_tips",
            "aligned_length",
            "ungapped_nonempty_columns",
        ],
    )
    print(f"Wrote {args.out_fasta}")
    print(f"Wrote {args.out_blocks}")
    print(f"Core alignment length: {next(iter(lengths))}")
    print(f"Used blocks: {sum(1 for row in block_rows if row['used_in_core'] == 'true')}")


if __name__ == "__main__":
    main()
