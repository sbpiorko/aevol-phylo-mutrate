from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
from collections import OrderedDict
from pathlib import Path
from typing import Iterable


MUTATION_RATE_KEYS = (
    "POINT_MUTATION_RATE",
    "SMALL_INSERTION_RATE",
    "SMALL_DELETION_RATE",
)


def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def parse_scalar(value: str):
    value = value.strip()
    if value.lower() in {"true", "false"}:
        return value.lower() == "true"
    try:
        if re.search(r"[.eE]", value):
            return float(value)
        return int(value)
    except ValueError:
        return value


def read_simple_yaml(path: Path) -> OrderedDict:
    """Read the small key/value YAML subset used by this scaffold."""
    data = OrderedDict()
    current_key = None

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue

        if not line.startswith((" ", "\t")):
            if ":" not in line:
                raise ValueError(f"Expected 'key: value' in {path}: {raw_line}")
            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip()
            if value == "":
                data[key] = OrderedDict()
                current_key = key
            else:
                data[key] = parse_scalar(value)
                current_key = key
            continue

        if current_key is None or not isinstance(data[current_key], OrderedDict):
            raise ValueError(f"Unexpected nested YAML line in {path}: {raw_line}")
        if ":" not in line:
            raise ValueError(f"Expected nested 'key: value' in {path}: {raw_line}")
        key, value = line.split(":", 1)
        data[current_key][key.strip()] = parse_scalar(value.strip())

    return data


def read_run_config(root: Path | None = None) -> OrderedDict:
    root = root or project_root()
    return read_simple_yaml(root / "config" / "run_config.yaml")


def read_mutation_grid(root: Path | None = None) -> OrderedDict[str, float]:
    root = root or project_root()
    raw = read_simple_yaml(root / "config" / "mutation_grid.yaml")
    grid: OrderedDict[str, float] = OrderedDict()
    for condition, value in raw.items():
        if isinstance(value, dict):
            multiplier = value.get("multiplier")
        else:
            multiplier = value
        if multiplier is None:
            raise ValueError(f"Missing multiplier for mutation condition {condition}")
        grid[str(condition)] = float(multiplier)
    return grid


def rel_path(root: Path, maybe_relative: str | Path) -> Path:
    path = Path(maybe_relative)
    if path.is_absolute():
        return path
    return root / path


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def comma_list(value) -> list[str]:
    if isinstance(value, list):
        return [str(v).strip() for v in value if str(v).strip()]
    return [part.strip() for part in str(value).split(",") if part.strip()]


def resolve_executable(candidates: Iterable[str]) -> str:
    for name in candidates:
        found = shutil.which(str(name))
        if found:
            return found
    raise FileNotFoundError(
        "None of these executables were found on PATH: "
        + ", ".join(str(c) for c in candidates)
    )


def run_logged(args: list[str], cwd: Path, log_path: Path) -> None:
    ensure_dir(log_path.parent)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("$ " + " ".join(args) + "\n\n")
        log.flush()
        proc = subprocess.run(
            args,
            cwd=str(cwd),
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        log.write(f"\n[exit_code] {proc.returncode}\n")
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {proc.returncode}. See {log_path}"
        )


def split_comment(line: str) -> tuple[str, str]:
    if "#" not in line:
        return line.rstrip("\n"), ""
    before, comment = line.rstrip("\n").split("#", 1)
    return before.rstrip(), " #" + comment


def first_numeric_token(tokens: list[str]) -> int | None:
    for i, token in enumerate(tokens):
        try:
            float(token)
            return i
        except ValueError:
            continue
    return None


def patch_param_text(
    text: str,
    replacements: dict[str, float | int | str],
    *,
    append_missing: bool = False,
) -> tuple[str, dict[str, tuple[str | None, str]]]:
    """Patch Aevol parameter lines while preserving comments."""
    replacements = {key.upper(): value for key, value in replacements.items()}
    seen: set[str] = set()
    changes: dict[str, tuple[str | None, str]] = {}
    out_lines: list[str] = []

    for raw_line in text.splitlines():
        body, comment = split_comment(raw_line)
        stripped = body.strip()
        if not stripped:
            out_lines.append(raw_line)
            continue

        tokens = stripped.split()
        key = tokens[0].upper()
        if key not in replacements:
            out_lines.append(raw_line)
            continue

        value = str(replacements[key])
        idx = first_numeric_token(tokens[1:])
        if idx is None:
            if len(tokens) == 1:
                old = None
                tokens.append(value)
            else:
                old = tokens[1]
                tokens[1] = value
        else:
            actual_idx = idx + 1
            old = tokens[actual_idx]
            tokens[actual_idx] = value

        seen.add(key)
        changes[key] = (old, value)
        indent = raw_line[: len(raw_line) - len(raw_line.lstrip())]
        out_lines.append(indent + " ".join(tokens) + comment)

    missing = [key for key in replacements if key not in seen]
    if missing and not append_missing:
        raise KeyError("Missing parameter(s): " + ", ".join(missing))

    for key in missing:
        value = str(replacements[key])
        out_lines.append(f"{key} {value}")
        changes[key] = (None, value)

    return "\n".join(out_lines) + "\n", changes


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks).upper()))
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line)

    if header is not None:
        records.append((header, "".join(chunks).upper()))

    return records


def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 80) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for i in range(0, len(sequence), width):
                handle.write(sequence[i : i + width] + "\n")


def first_fasta_record(path: Path) -> tuple[str, str]:
    records = parse_fasta(path)
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records[0]


def find_candidate_fastas(directory: Path) -> list[Path]:
    patterns = ("*.fa", "*.fasta", "*.fas", "*.fna")
    files: list[Path] = []
    for pattern in patterns:
        files.extend(directory.glob(pattern))
    return sorted(set(files), key=lambda p: p.stat().st_mtime, reverse=True)


def extract_best_genome(
    *,
    root: Path,
    run_dir: Path,
    generation: int,
    output_fasta: Path,
    config: OrderedDict,
    log_name: str = "extract.log",
) -> Path:
    extract_cmd = resolve_executable(comma_list(config["aevol_extract_candidates"]))
    cmd_name = Path(extract_cmd).name.lower()
    output_fasta_abs = output_fasta.resolve()

    if "misc_extract" in cmd_name:
        args = [extract_cmd, "-t", str(generation), "-S", str(output_fasta_abs)]
        run_logged(args, run_dir, run_dir / log_name)
        if not output_fasta.exists():
            raise FileNotFoundError(f"Expected extracted FASTA at {output_fasta}")
        return output_fasta

    before = set(find_candidate_fastas(run_dir))
    args = [extract_cmd, "-b", str(generation), "-e", str(generation)]
    run_logged(args, run_dir, run_dir / log_name)
    after = set(find_candidate_fastas(run_dir))
    new_files = sorted(after - before, key=lambda p: p.stat().st_mtime, reverse=True)
    candidates = new_files or find_candidate_fastas(run_dir)

    for candidate in candidates:
        if candidate.resolve() == output_fasta_abs:
            return output_fasta
        if candidate.name in {"parent.fa", "params.in"}:
            continue
        records = parse_fasta(candidate)
        if records:
            shutil.copyfile(candidate, output_fasta)
            return output_fasta

    raise FileNotFoundError(
        f"Could not identify an extracted FASTA in {run_dir}. "
        "If your Aevol extractor needs different flags, adjust "
        "extract_best_genome() in scripts/pipeline_common.py."
    )


def patch_seed_in_file(path: Path, seed: int) -> None:
    text = path.read_text(encoding="utf-8")
    patched, _ = patch_param_text(text, {"SEED": seed}, append_missing=True)
    path.write_text(patched, encoding="utf-8")

