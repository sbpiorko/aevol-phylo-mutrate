from __future__ import annotations

from pathlib import Path

from pipeline_common import (
    MUTATION_RATE_KEYS,
    comma_list,
    project_root,
    read_mutation_grid,
    read_run_config,
    rel_path,
    resolve_executable,
)


def main() -> None:
    root = project_root()
    config = read_run_config(root)
    grid = read_mutation_grid(root)
    base_params = rel_path(root, config["base_params"])

    print(f"Project root: {root}")
    print(f"Mutation conditions: {', '.join(grid.keys())}")

    if not base_params.exists():
        raise FileNotFoundError(base_params)

    text = base_params.read_text(encoding="utf-8")
    missing_keys = [key for key in MUTATION_RATE_KEYS if key not in text]
    if missing_keys:
        raise ValueError(
            f"{base_params} is missing required rate keys: {', '.join(missing_keys)}"
        )
    print(f"Base params: {base_params}")

    for name, key in [
        ("create", "aevol_create_candidates"),
        ("run", "aevol_run_candidates"),
        ("extract", "aevol_extract_candidates"),
    ]:
        try:
            exe = resolve_executable(comma_list(config[key]))
            print(f"Aevol {name}: {exe}")
        except FileNotFoundError as exc:
            print(f"WARNING: {exc}")

    print("Setup check finished.")


if __name__ == "__main__":
    main()

