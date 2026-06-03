#!/usr/bin/env python3
"""Score inferred phylogenies against the known simulated tree.

This script intentionally implements topology-first scoring without external
tree libraries so it can run in the lightweight Codex environment as well as
inside WSL. Branch-length-sensitive metrics are reported only when branch
lengths pass basic validation.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, field
from pathlib import Path

from pipeline_common import ensure_dir, project_root, read_mutation_grid, write_tsv


@dataclass(eq=False)
class Node:
    name: str = ""
    length: float | None = None
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None

    def is_leaf(self) -> bool:
        return not self.children


class NewickParser:
    def __init__(self, text: str) -> None:
        self.text = text.strip()
        if self.text.endswith(";"):
            self.text = self.text[:-1]
        self.index = 0

    def parse(self) -> Node:
        node = self._parse_subtree()
        self._skip_ws()
        if self.index != len(self.text):
            raise ValueError(f"Unexpected Newick content at position {self.index}")
        self._set_parents(node)
        return node

    def _parse_subtree(self) -> Node:
        self._skip_ws()
        if self._peek() == "(":
            self.index += 1
            children = []
            while True:
                children.append(self._parse_subtree())
                self._skip_ws()
                char = self._peek()
                if char == ",":
                    self.index += 1
                    continue
                if char == ")":
                    self.index += 1
                    break
                raise ValueError(f"Expected ',' or ')' at position {self.index}")
            name = self._parse_label()
            length = self._parse_length()
            return Node(name=name, length=length, children=children)

        name = self._parse_label()
        if not name:
            raise ValueError(f"Expected leaf label at position {self.index}")
        length = self._parse_length()
        return Node(name=name, length=length)

    def _parse_label(self) -> str:
        self._skip_ws()
        start = self.index
        while self.index < len(self.text) and self.text[self.index] not in ":,()":
            self.index += 1
        return self.text[start : self.index].strip()

    def _parse_length(self) -> float | None:
        self._skip_ws()
        if self._peek() != ":":
            return None
        self.index += 1
        self._skip_ws()
        start = self.index
        while self.index < len(self.text) and self.text[self.index] not in ",()":
            self.index += 1
        raw = self.text[start : self.index].strip()
        if not raw:
            raise ValueError(f"Missing branch length at position {start}")
        return float(raw)

    def _peek(self) -> str:
        if self.index >= len(self.text):
            return ""
        return self.text[self.index]

    def _skip_ws(self) -> None:
        while self.index < len(self.text) and self.text[self.index].isspace():
            self.index += 1

    def _set_parents(self, node: Node) -> None:
        for child in node.children:
            child.parent = node
            self._set_parents(child)


def parse_newick(path: Path) -> Node:
    return NewickParser(path.read_text(encoding="utf-8")).parse()


def leaves(node: Node) -> list[Node]:
    if node.is_leaf():
        return [node]
    out: list[Node] = []
    for child in node.children:
        out.extend(leaves(child))
    return out


def leaf_names(node: Node) -> set[str]:
    return {leaf.name for leaf in leaves(node)}


def collect_lengths(node: Node) -> list[float]:
    lengths = []
    for child in node.children:
        if child.length is not None:
            lengths.append(child.length)
        lengths.extend(collect_lengths(child))
    return lengths


def descendant_tips(node: Node) -> set[str]:
    if node.is_leaf():
        return {node.name}
    tips: set[str] = set()
    for child in node.children:
        tips.update(descendant_tips(child))
    return tips


def canonical_split(side: set[str], all_tips: set[str]) -> tuple[str, ...]:
    other = all_tips - side
    chosen = side
    if len(other) < len(side) or (len(other) == len(side) and sorted(other) < sorted(side)):
        chosen = other
    return tuple(sorted(chosen))


def unrooted_splits(root: Node) -> set[tuple[str, ...]]:
    all_tips = leaf_names(root)
    splits: set[tuple[str, ...]] = set()

    def walk(node: Node) -> None:
        for child in node.children:
            side = descendant_tips(child)
            other = all_tips - side
            if len(side) > 1 and len(other) > 1:
                splits.add(canonical_split(side, all_tips))
            walk(child)

    walk(root)
    return splits


def format_split(split: tuple[str, ...]) -> str:
    return ",".join(split)


def tip_distance_vectors(root: Node, ordered_tips: list[str]) -> dict[str, dict[str, float]]:
    paths: dict[str, dict[Node, float]] = {}

    def walk(node: Node, lineage: dict[Node, float], distance: float) -> None:
        next_lineage = dict(lineage)
        next_lineage[node] = distance
        if node.is_leaf():
            paths[node.name] = next_lineage
            return
        for child in node.children:
            walk(child, next_lineage, distance + float(child.length or 0.0))

    walk(root, {}, 0.0)

    distances: dict[str, dict[str, float]] = {tip: {} for tip in ordered_tips}
    for i, left in enumerate(ordered_tips):
        for right in ordered_tips[i + 1 :]:
            common = set(paths[left]).intersection(paths[right])
            lca_distance = max(paths[left][node] for node in common)
            distance = paths[left][root]  # quiet type-checker; overwritten below
            distance = paths[left][next(iter(paths[left]))] if False else (
                max(paths[left].values()) + max(paths[right].values()) - 2 * lca_distance
            )
            distances[left][right] = distance
            distances[right][left] = distance
    return distances


def pearson(left: list[float], right: list[float]) -> float | None:
    if len(left) != len(right) or len(left) < 2:
        return None
    mean_left = sum(left) / len(left)
    mean_right = sum(right) / len(right)
    numerator = sum((a - mean_left) * (b - mean_right) for a, b in zip(left, right))
    denom_left = math.sqrt(sum((a - mean_left) ** 2 for a in left))
    denom_right = math.sqrt(sum((b - mean_right) ** 2 for b in right))
    denominator = denom_left * denom_right
    if denominator == 0:
        return None
    return numerator / denominator


def patristic_correlation(true_root: Node, inferred_root: Node, ordered_tips: list[str]) -> float | None:
    true_dist = tip_distance_vectors(true_root, ordered_tips)
    inferred_dist = tip_distance_vectors(inferred_root, ordered_tips)
    true_values: list[float] = []
    inferred_values: list[float] = []
    for i, left in enumerate(ordered_tips):
        for right in ordered_tips[i + 1 :]:
            true_values.append(true_dist[left][right])
            inferred_values.append(inferred_dist[left][right])
    return pearson(true_values, inferred_values)


def parse_args() -> argparse.Namespace:
    root = project_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--true-tree", type=Path, default=root / "data" / "true_trees" / "true_tree_16tips.nwk")
    parser.add_argument(
        "--methods",
        nargs="+",
        default=["mash_bionj", "mafft_iqtree", "mauve_iqtree"],
        help="Phylogeny method result folders to score.",
    )
    parser.add_argument("--conditions", nargs="+", default=list(read_mutation_grid(root).keys()))
    parser.add_argument("--replicates", nargs="+", default=["rep_001"])
    return parser.parse_args()


def main() -> None:
    root = project_root()
    args = parse_args()
    true_tree_path = args.true_tree if args.true_tree.is_absolute() else root / args.true_tree
    true_root = parse_newick(true_tree_path)
    true_tips = leaf_names(true_root)
    ordered_tips = sorted(true_tips)
    true_splits = unrooted_splits(true_root)
    max_rf = 2 * (len(true_tips) - 3) if len(true_tips) > 3 else 0

    score_rows: list[dict] = []
    validation_rows: list[dict] = []
    split_rows: list[dict] = []
    patristic_rows: list[dict] = []

    for method in args.methods:
        for condition in args.conditions:
            for replicate in args.replicates:
                tree_path = (
                    root
                    / "results"
                    / "phylo"
                    / method
                    / condition
                    / replicate
                    / "inferred_tree.nwk"
                )
                row_base = {
                    "method": method,
                    "condition": condition,
                    "replicate": replicate,
                    "tree": tree_path.relative_to(root).as_posix(),
                }

                if not tree_path.exists():
                    validation_rows.append(
                        {
                            **row_base,
                            "status": "missing_tree",
                            "tip_count": 0,
                            "missing_tips": "",
                            "extra_tips": "",
                            "branch_count": 0,
                            "negative_branch_lengths": 0,
                            "min_branch_length": "",
                            "max_branch_length": "",
                        }
                    )
                    continue

                inferred_root = parse_newick(tree_path)
                inferred_tips = leaf_names(inferred_root)
                missing_tips = sorted(true_tips - inferred_tips)
                extra_tips = sorted(inferred_tips - true_tips)
                lengths = collect_lengths(inferred_root)
                negative_lengths = [length for length in lengths if length < -1e-12]
                status = "ok"
                if missing_tips or extra_tips:
                    status = "tip_mismatch"
                elif negative_lengths:
                    status = "ok_with_negative_branch_lengths"

                validation_rows.append(
                    {
                        **row_base,
                        "status": status,
                        "tip_count": len(inferred_tips),
                        "missing_tips": ",".join(missing_tips),
                        "extra_tips": ",".join(extra_tips),
                        "branch_count": len(lengths),
                        "negative_branch_lengths": len(negative_lengths),
                        "min_branch_length": f"{min(lengths):.12g}" if lengths else "",
                        "max_branch_length": f"{max(lengths):.12g}" if lengths else "",
                    }
                )

                if missing_tips or extra_tips:
                    continue

                inferred_splits = unrooted_splits(inferred_root)
                matching = true_splits & inferred_splits
                missing = true_splits - inferred_splits
                extra = inferred_splits - true_splits
                rf_distance = len(missing) + len(extra)
                normalized_rf = rf_distance / max_rf if max_rf else 0.0

                score_rows.append(
                    {
                        **row_base,
                        "true_split_count": len(true_splits),
                        "inferred_split_count": len(inferred_splits),
                        "matching_splits": len(matching),
                        "missing_splits": len(missing),
                        "extra_splits": len(extra),
                        "rf_distance": rf_distance,
                        "normalized_rf": f"{normalized_rf:.6f}",
                        "negative_branch_lengths": len(negative_lengths),
                    }
                )

                for split_type, splits in [
                    ("matching", matching),
                    ("missing_from_inferred", missing),
                    ("extra_in_inferred", extra),
                ]:
                    for split in sorted(splits):
                        split_rows.append({**row_base, "split_type": split_type, "split": format_split(split)})

                if negative_lengths:
                    patristic_rows.append(
                        {
                            **row_base,
                            "status": "skipped_negative_branch_lengths",
                            "patristic_correlation": "",
                        }
                    )
                else:
                    corr = patristic_correlation(true_root, inferred_root, ordered_tips)
                    patristic_rows.append(
                        {
                            **row_base,
                            "status": "ok" if corr is not None else "undefined",
                            "patristic_correlation": f"{corr:.6f}" if corr is not None else "",
                        }
                    )

                print(
                    f"{method} {condition} {replicate}: RF={rf_distance}, "
                    f"normalized RF={normalized_rf:.3f}, status={status}"
                )

    out_dir = ensure_dir(root / "results" / "tree_scores")
    write_tsv(
        out_dir / "tree_scores.tsv",
        score_rows,
        [
            "method",
            "condition",
            "replicate",
            "tree",
            "true_split_count",
            "inferred_split_count",
            "matching_splits",
            "missing_splits",
            "extra_splits",
            "rf_distance",
            "normalized_rf",
            "negative_branch_lengths",
        ],
    )
    write_tsv(
        out_dir / "tree_validation.tsv",
        validation_rows,
        [
            "method",
            "condition",
            "replicate",
            "tree",
            "status",
            "tip_count",
            "missing_tips",
            "extra_tips",
            "branch_count",
            "negative_branch_lengths",
            "min_branch_length",
            "max_branch_length",
        ],
    )
    write_tsv(
        out_dir / "tree_splits.tsv",
        split_rows,
        ["method", "condition", "replicate", "tree", "split_type", "split"],
    )
    write_tsv(
        out_dir / "patristic_correlations.tsv",
        patristic_rows,
        ["method", "condition", "replicate", "tree", "status", "patristic_correlation"],
    )

    print(f"Wrote {out_dir / 'tree_scores.tsv'}")
    print(f"Wrote {out_dir / 'tree_validation.tsv'}")
    print(f"Wrote {out_dir / 'tree_splits.tsv'}")
    print(f"Wrote {out_dir / 'patristic_correlations.tsv'}")


if __name__ == "__main__":
    main()
