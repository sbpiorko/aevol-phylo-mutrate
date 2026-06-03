from __future__ import annotations

import argparse
import shutil

from pipeline_common import (
    comma_list,
    ensure_dir,
    extract_best_genome,
    patch_seed_in_file,
    project_root,
    read_run_config,
    rel_path,
    resolve_executable,
    run_logged,
)


def parse_args() -> argparse.Namespace:
    root = project_root()
    config = read_run_config(root)
    parser = argparse.ArgumentParser(description="Pre-evolve one baseline ancestor.")
    parser.add_argument("--generations", type=int, default=int(config["ancestor_generations"]))
    parser.add_argument("--threads", type=int, default=int(config["threads"]))
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    root = project_root()
    config = read_run_config(root)
    args = parse_args()

    run_dir = ensure_dir(root / "runs" / "ancestor")
    out_fasta = root / "data" / "ancestor" / "ancestor_best.fa"
    ensure_dir(out_fasta.parent)

    if out_fasta.exists() and not args.force:
        print(f"Ancestor FASTA already exists: {out_fasta}")
        print("Use --force to rerun ancestor pre-evolution.")
        return

    base_params = rel_path(root, config["base_params"])
    params = run_dir / "params.in"
    shutil.copyfile(base_params, params)
    patch_seed_in_file(params, int(config["root_seed"]))

    create_cmd = resolve_executable(comma_list(config["aevol_create_candidates"]))
    run_cmd = resolve_executable(comma_list(config["aevol_run_candidates"]))

    run_logged([create_cmd, "params.in"], run_dir, run_dir / "create.log")
    run_logged(
        [run_cmd, "-e", str(args.generations), "-p", str(args.threads)],
        run_dir,
        run_dir / "run.log",
    )
    extract_best_genome(
        root=root,
        run_dir=run_dir,
        generation=args.generations,
        output_fasta=out_fasta,
        config=config,
    )

    metadata = root / "data" / "ancestor" / "ancestor_metadata.tsv"
    metadata.write_text(
        "key\tvalue\n"
        f"generations\t{args.generations}\n"
        f"seed\t{config['root_seed']}\n"
        f"params\t{params.relative_to(root).as_posix()}\n"
        f"fasta\t{out_fasta.relative_to(root).as_posix()}\n",
        encoding="utf-8",
    )

    print(f"Wrote {out_fasta}")
    print(f"Wrote {metadata}")


if __name__ == "__main__":
    main()

