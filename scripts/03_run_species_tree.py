from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from pipeline_common import (
    comma_list,
    ensure_dir,
    extract_best_genome,
    patch_seed_in_file,
    project_root,
    read_mutation_grid,
    read_run_config,
    read_tsv,
    rel_path,
    resolve_executable,
    run_logged,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    root = project_root()
    config = read_run_config(root)
    parser = argparse.ArgumentParser(
        description="Run the known species tree for one or more mutation conditions."
    )
    parser.add_argument(
        "--conditions",
        nargs="+",
        default=None,
        help="Subset of conditions to run. Default: all conditions in mutation_grid.yaml.",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=int(config["replicates"]),
        help="Number of replicate trees to run from rep_001 upward.",
    )
    parser.add_argument("--threads", type=int, default=int(config["threads"]))
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def seed_for(condition_index: int, replicate: int, edge_number: int, seed_base: int) -> int:
    return seed_base + (condition_index + 1) * 10000 + replicate * 1000 + edge_number


def branch_dir_name(edge: dict[str, str]) -> str:
    return f"{edge['edge_id']}_{edge['parent']}_to_{edge['child']}"


def run_branch(
    *,
    root: Path,
    config: dict,
    condition: str,
    condition_index: int,
    replicate: int,
    edge: dict[str, str],
    params_template: Path,
    parent_fasta: Path,
    create_cmd: str,
    run_cmd: str,
    threads: int,
    force: bool,
) -> tuple[dict, dict]:
    rep_label = f"rep_{replicate:03d}"
    run_dir = ensure_dir(root / "runs" / condition / rep_label / branch_dir_name(edge))
    child = edge["child"]
    child_fasta = root / "data" / "node_fastas" / condition / rep_label / f"{child}.fa"
    branch_child_fasta = run_dir / "final_child.fa"

    edge_number = int(edge["edge_id"].lstrip("E"))
    seed = seed_for(condition_index, replicate, edge_number, int(config["seed_base"]))

    seed_row = {
        "condition": condition,
        "replicate": replicate,
        "edge_id": edge["edge_id"],
        "parent": edge["parent"],
        "child": child,
        "seed": seed,
    }

    if child_fasta.exists() and branch_child_fasta.exists() and not force:
        status_row = {
            **seed_row,
            "run_dir": run_dir.relative_to(root).as_posix(),
            "status": "skipped_existing",
            "child_fasta": child_fasta.relative_to(root).as_posix(),
        }
        return seed_row, status_row

    if not parent_fasta.exists():
        raise FileNotFoundError(f"Missing parent FASTA for {edge['edge_id']}: {parent_fasta}")

    shutil.copyfile(params_template, run_dir / "params.in")
    patch_seed_in_file(run_dir / "params.in", seed)
    shutil.copyfile(parent_fasta, run_dir / "parent.fa")

    run_logged(
        [create_cmd, "params.in", "--fasta", "parent.fa"],
        run_dir,
        run_dir / "create.log",
    )
    generations = int(edge["branch_generations"])
    run_logged(
        [run_cmd, "-e", str(generations), "-p", str(threads)],
        run_dir,
        run_dir / "run.log",
    )

    extract_best_genome(
        root=root,
        run_dir=run_dir,
        generation=generations,
        output_fasta=branch_child_fasta,
        config=config,
    )
    ensure_dir(child_fasta.parent)
    shutil.copyfile(branch_child_fasta, child_fasta)

    status_row = {
        **seed_row,
        "run_dir": run_dir.relative_to(root).as_posix(),
        "status": "completed",
        "child_fasta": child_fasta.relative_to(root).as_posix(),
    }
    return seed_row, status_row


def merge_manifest_rows(
    existing_path: Path,
    new_rows: list[dict],
    key_fields: list[str],
) -> list[dict]:
    merged: dict[tuple[str, ...], dict] = {}
    if existing_path.exists():
        for row in read_tsv(existing_path):
            key = tuple(str(row[field]) for field in key_fields)
            merged[key] = row
    for row in new_rows:
        key = tuple(str(row[field]) for field in key_fields)
        merged[key] = row
    return list(merged.values())


def main() -> None:
    root = project_root()
    config = read_run_config(root)
    grid = read_mutation_grid(root)
    args = parse_args()

    conditions = args.conditions or list(grid.keys())
    unknown = [condition for condition in conditions if condition not in grid]
    if unknown:
        raise ValueError(f"Unknown condition(s): {', '.join(unknown)}")

    edges_path = rel_path(root, config["tree_edges"])
    edges = read_tsv(edges_path)
    if not edges:
        raise ValueError(f"No edges found in {edges_path}. Run 00_make_tree.py first.")

    ancestor_fasta = root / "data" / "ancestor" / "ancestor_best.fa"
    if not ancestor_fasta.exists():
        raise FileNotFoundError(
            f"Missing ancestor FASTA: {ancestor_fasta}. Run 02_pre_evolve_ancestor.py first."
        )

    create_cmd = resolve_executable(comma_list(config["aevol_create_candidates"]))
    run_cmd = resolve_executable(comma_list(config["aevol_run_candidates"]))

    seed_rows: list[dict] = []
    status_rows: list[dict] = []

    for condition in conditions:
        condition_index = list(grid.keys()).index(condition)
        params_template = (
            rel_path(root, config["generated_params_dir"]) / f"{condition}_params.in"
        )
        if not params_template.exists():
            raise FileNotFoundError(
                f"Missing generated params for {condition}: {params_template}. "
                "Run 01_make_mutation_params.py first."
            )

        for replicate in range(1, args.replicates + 1):
            rep_label = f"rep_{replicate:03d}"
            node_dir = ensure_dir(root / "data" / "node_fastas" / condition / rep_label)
            root_fasta = node_dir / "root.fa"
            if args.force or not root_fasta.exists():
                shutil.copyfile(ancestor_fasta, root_fasta)

            node_fastas = {"root": root_fasta}
            for edge in edges:
                parent_fasta = node_fastas.get(edge["parent"])
                if parent_fasta is None:
                    parent_fasta = node_dir / f"{edge['parent']}.fa"
                seed_row, status_row = run_branch(
                    root=root,
                    config=config,
                    condition=condition,
                    condition_index=condition_index,
                    replicate=replicate,
                    edge=edge,
                    params_template=params_template,
                    parent_fasta=parent_fasta,
                    create_cmd=create_cmd,
                    run_cmd=run_cmd,
                    threads=args.threads,
                    force=args.force,
                )
                seed_rows.append(seed_row)
                status_rows.append(status_row)
                node_fastas[edge["child"]] = (
                    root
                    / "data"
                    / "node_fastas"
                    / condition
                    / rep_label
                    / f"{edge['child']}.fa"
                )
                print(
                    f"{condition} {rep_label} {edge['edge_id']} "
                    f"{edge['parent']}->{edge['child']}: {status_row['status']}"
                )

    seed_manifest = root / "data" / "manifests" / "branch_seeds.tsv"
    status_manifest = root / "data" / "manifests" / "simulation_status.tsv"
    seed_rows = merge_manifest_rows(
        seed_manifest,
        seed_rows,
        ["condition", "replicate", "edge_id"],
    )
    status_rows = merge_manifest_rows(
        status_manifest,
        status_rows,
        ["condition", "replicate", "edge_id"],
    )

    write_tsv(
        seed_manifest,
        seed_rows,
        ["condition", "replicate", "edge_id", "parent", "child", "seed"],
    )
    write_tsv(
        status_manifest,
        status_rows,
        [
            "condition",
            "replicate",
            "edge_id",
            "parent",
            "child",
            "seed",
            "run_dir",
            "status",
            "child_fasta",
        ],
    )

    print("Species-tree simulations finished.")


if __name__ == "__main__":
    main()
