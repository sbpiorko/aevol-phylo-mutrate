from __future__ import annotations

import argparse
from pathlib import Path

from pipeline_common import (
    MUTATION_RATE_KEYS,
    ensure_dir,
    patch_param_text,
    project_root,
    read_mutation_grid,
    read_run_config,
    rel_path,
    write_tsv,
)


def numeric_value_for_key(text: str, key: str) -> float:
    for raw_line in text.splitlines():
        body = raw_line.split("#", 1)[0].strip()
        if not body:
            continue
        tokens = body.split()
        if tokens and tokens[0].upper() == key.upper():
            for token in tokens[1:]:
                try:
                    return float(token)
                except ValueError:
                    continue
    raise KeyError(f"Could not find numeric value for {key}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate condition-specific Aevol params by scaling local rates."
    )
    parser.add_argument("--base-params", default=None)
    parser.add_argument("--out-dir", default=None)
    return parser.parse_args()


def main() -> None:
    root = project_root()
    config = read_run_config(root)
    args = parse_args()

    base_params = rel_path(root, args.base_params or config["base_params"])
    out_dir = rel_path(root, args.out_dir or config["generated_params_dir"])
    ensure_dir(out_dir)

    base_text = base_params.read_text(encoding="utf-8")
    base_values = {key: numeric_value_for_key(base_text, key) for key in MUTATION_RATE_KEYS}
    grid = read_mutation_grid(root)

    manifest_rows = []
    for condition, multiplier in grid.items():
        replacements = {
            key: f"{base_values[key] * multiplier:.12g}" for key in MUTATION_RATE_KEYS
        }
        patched, changes = patch_param_text(base_text, replacements)
        out_path = out_dir / f"{condition}_params.in"
        out_path.write_text(patched, encoding="utf-8")

        for key in MUTATION_RATE_KEYS:
            old, new = changes[key]
            manifest_rows.append(
                {
                    "condition": condition,
                    "multiplier": multiplier,
                    "parameter": key,
                    "base_value": old,
                    "scaled_value": new,
                    "params_file": out_path.relative_to(root).as_posix(),
                }
            )

    manifest = root / "data" / "manifests" / "mutation_params_resolved.tsv"
    write_tsv(
        manifest,
        manifest_rows,
        [
            "condition",
            "multiplier",
            "parameter",
            "base_value",
            "scaled_value",
            "params_file",
        ],
    )

    print(f"Wrote generated params to {out_dir}")
    print(f"Wrote {manifest}")


if __name__ == "__main__":
    main()

