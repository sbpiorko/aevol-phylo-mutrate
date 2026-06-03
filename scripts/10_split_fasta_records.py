#!/usr/bin/env python3
"""Split a multi-FASTA into one FASTA file per record."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from pipeline_common import ensure_dir, parse_fasta, write_fasta


SAFE_NAME = re.compile(r"^[A-Za-z0-9_.-]+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    records = parse_fasta(args.input)
    if not records:
        raise ValueError(f"No FASTA records found in {args.input}")

    out_dir = ensure_dir(args.out_dir)
    written = []
    for header, sequence in records:
        record_id = header.split()[0]
        if not SAFE_NAME.match(record_id):
            raise ValueError(f"Unsafe FASTA record id for filename: {record_id!r}")
        out_path = out_dir / f"{record_id}.fa"
        write_fasta(out_path, [(record_id, sequence)])
        written.append(out_path)

    manifest = out_dir / "split_fastas.txt"
    manifest.write_text(
        "\n".join(path.as_posix() for path in written) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {len(written)} split FASTAs to {out_dir}")
    print(f"Wrote {manifest}")


if __name__ == "__main__":
    main()
