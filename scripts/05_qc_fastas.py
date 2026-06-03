from __future__ import annotations

import argparse
from collections import Counter

from pipeline_common import (
    parse_fasta,
    project_root,
    read_mutation_grid,
    read_run_config,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    root = project_root()
    config = read_run_config(root)
    parser = argparse.ArgumentParser(description="QC final tip FASTA files.")
    parser.add_argument("--conditions", nargs="+", default=None)
    parser.add_argument("--replicates", type=int, default=int(config["replicates"]))
    return parser.parse_args()


def main() -> None:
    root = project_root()
    config = read_run_config(root)
    grid = read_mutation_grid(root)
    args = parse_args()

    expected_tips = int(config["tips"])
    expected_names = {f"S{i:02d}" for i in range(1, expected_tips + 1)}
    allowed = set("ATGCN")

    fasta_qc_rows = []
    length_rows = []
    status_rows = []

    for condition in args.conditions or list(grid.keys()):
        for replicate in range(1, args.replicates + 1):
            rep_label = f"rep_{replicate:03d}"
            all_tips = (
                root
                / "data"
                / "final_fastas"
                / condition
                / rep_label
                / "all_tips.fa"
            )
            status = "ok"
            message = ""

            if not all_tips.exists():
                status_rows.append(
                    {
                        "condition": condition,
                        "replicate": replicate,
                        "status": "missing",
                        "message": f"Missing {all_tips}",
                    }
                )
                continue

            records = parse_fasta(all_tips)
            names = [header.split()[0] for header, _ in records]
            sequences = [sequence for _, sequence in records]
            name_counts = Counter(names)
            duplicate_names = sorted(name for name, count in name_counts.items() if count > 1)
            duplicate_sequences = len(sequences) - len(set(sequences))
            missing_names = sorted(expected_names - set(names))
            extra_names = sorted(set(names) - expected_names)
            invalid_records = [
                name for name, sequence in zip(names, sequences) if set(sequence) - allowed
            ]
            lengths = [len(sequence) for sequence in sequences]

            if len(records) != expected_tips:
                status = "fail"
                message += f"expected {expected_tips} records, found {len(records)}; "
            if missing_names:
                status = "fail"
                message += "missing " + ",".join(missing_names) + "; "
            if extra_names:
                status = "fail"
                message += "extra " + ",".join(extra_names) + "; "
            if duplicate_names:
                status = "fail"
                message += "duplicate names " + ",".join(duplicate_names) + "; "
            if invalid_records:
                status = "fail"
                message += "invalid bases in " + ",".join(invalid_records) + "; "

            for name, sequence in zip(names, sequences):
                length_rows.append(
                    {
                        "condition": condition,
                        "replicate": replicate,
                        "tip": name,
                        "length": len(sequence),
                    }
                )

            fasta_qc_rows.append(
                {
                    "condition": condition,
                    "replicate": replicate,
                    "records": len(records),
                    "min_length": min(lengths) if lengths else 0,
                    "max_length": max(lengths) if lengths else 0,
                    "mean_length": f"{sum(lengths) / len(lengths):.2f}" if lengths else "0",
                    "duplicate_sequences": duplicate_sequences,
                    "invalid_records": len(invalid_records),
                    "fasta": all_tips.relative_to(root).as_posix(),
                }
            )
            status_rows.append(
                {
                    "condition": condition,
                    "replicate": replicate,
                    "status": status,
                    "message": message.strip(),
                }
            )

    out_dir = root / "results" / "simulation_qc"
    write_tsv(
        out_dir / "fasta_qc.tsv",
        fasta_qc_rows,
        [
            "condition",
            "replicate",
            "records",
            "min_length",
            "max_length",
            "mean_length",
            "duplicate_sequences",
            "invalid_records",
            "fasta",
        ],
    )
    write_tsv(
        out_dir / "genome_lengths.tsv",
        length_rows,
        ["condition", "replicate", "tip", "length"],
    )
    write_tsv(
        out_dir / "simulation_status.tsv",
        status_rows,
        ["condition", "replicate", "status", "message"],
    )

    failed = [row for row in status_rows if row["status"] != "ok"]
    if failed:
        print(f"QC finished with {len(failed)} failing condition/replicate combinations.")
        raise SystemExit(1)

    print("QC passed for all final FASTA files.")


if __name__ == "__main__":
    main()

