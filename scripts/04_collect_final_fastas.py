from __future__ import annotations

import argparse

from pipeline_common import (
    ensure_dir,
    first_fasta_record,
    project_root,
    read_mutation_grid,
    read_run_config,
    read_tsv,
    rel_path,
    write_fasta,
)


def parse_args() -> argparse.Namespace:
    root = project_root()
    config = read_run_config(root)
    parser = argparse.ArgumentParser(description="Collect final tip FASTAs.")
    parser.add_argument("--conditions", nargs="+", default=None)
    parser.add_argument("--replicates", type=int, default=int(config["replicates"]))
    return parser.parse_args()


def main() -> None:
    root = project_root()
    config = read_run_config(root)
    grid = read_mutation_grid(root)
    args = parse_args()

    conditions = args.conditions or list(grid.keys())
    edges = read_tsv(rel_path(root, config["tree_edges"]))
    tips = sorted(edge["child"] for edge in edges if edge["is_tip"] == "true")

    for condition in conditions:
        if condition not in grid:
            raise ValueError(f"Unknown condition: {condition}")
        for replicate in range(1, args.replicates + 1):
            rep_label = f"rep_{replicate:03d}"
            node_dir = root / "data" / "node_fastas" / condition / rep_label
            out_dir = ensure_dir(root / "data" / "final_fastas" / condition / rep_label)
            all_records = []

            for tip in tips:
                source = node_dir / f"{tip}.fa"
                if not source.exists():
                    raise FileNotFoundError(f"Missing tip FASTA: {source}")
                _, sequence = first_fasta_record(source)
                record = (tip, sequence)
                write_fasta(out_dir / f"{tip}.fa", [record])
                all_records.append(record)

            all_tips = out_dir / "all_tips.fa"
            write_fasta(all_tips, all_records)
            print(f"Wrote {all_tips} with {len(all_records)} tip genomes")


if __name__ == "__main__":
    main()

