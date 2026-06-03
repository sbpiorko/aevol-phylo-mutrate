#!/usr/bin/env python3
"""Create lightweight SVG summary plots for current phylogeny results."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from pipeline_common import ensure_dir, project_root


CONDITIONS = ["very_low", "low", "normal", "high", "very_high"]
CONDITION_LABELS = {
    "very_low": "Very low",
    "low": "Low",
    "normal": "Normal",
    "high": "High",
    "very_high": "Very high",
}
METHOD_LABELS = {
    "mash_bionj": "Mash/BioNJ",
    "mafft_iqtree": "MAFFT/IQ-TREE",
    "mauve_iqtree": "Mauve/IQ-TREE",
}
METHOD_COLORS = {
    "mash_bionj": "#2f6f73",
    "mafft_iqtree": "#b24f2a",
    "mauve_iqtree": "#5f5aa2",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def svg_text(x: float, y: float, text: str, *, size: int = 12, anchor: str = "middle", weight: str = "400") -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
        f'font-size="{size}" font-weight="{weight}" fill="#1d2528">{escape(text)}</text>'
    )


def escape(text: str) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def bar_chart(
    *,
    path: Path,
    title: str,
    rows: list[dict[str, str]],
    value_key: str,
    y_label: str,
    y_max: float | None = None,
) -> None:
    methods = sorted({row["method"] for row in rows}, key=lambda m: list(METHOD_LABELS).index(m) if m in METHOD_LABELS else 99)
    values = {(row["condition"], row["method"]): float(row[value_key]) for row in rows if row.get(value_key)}
    if y_max is None:
        y_max = max(values.values(), default=1.0)
        if y_max <= 0:
            y_max = 1.0

    width = 980
    height = 560
    left = 88
    right = 28
    top = 64
    bottom = 96
    plot_w = width - left - right
    plot_h = height - top - bottom
    group_w = plot_w / len(CONDITIONS)
    bar_gap = 8
    bar_w = min(42, (group_w - 34) / max(len(methods), 1) - bar_gap)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf7"/>',
        svg_text(width / 2, 34, title, size=20, weight="700"),
        svg_text(20, top + plot_h / 2, y_label, size=13, anchor="middle"),
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#2d3336" stroke-width="1.2"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#2d3336" stroke-width="1.2"/>',
    ]

    ticks = 5
    for i in range(ticks + 1):
        value = y_max * i / ticks
        y = top + plot_h - (value / y_max) * plot_h
        parts.append(f'<line x1="{left - 5}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" stroke="#d9d4ca" stroke-width="1"/>')
        parts.append(svg_text(left - 10, y + 4, f"{value:.2f}", size=11, anchor="end"))

    for c_index, condition in enumerate(CONDITIONS):
        group_x = left + c_index * group_w
        center = group_x + group_w / 2
        parts.append(svg_text(center, top + plot_h + 34, CONDITION_LABELS[condition], size=12))
        start_x = center - ((bar_w + bar_gap) * len(methods) - bar_gap) / 2
        for m_index, method in enumerate(methods):
            value = values.get((condition, method))
            if value is None:
                continue
            bar_h = (value / y_max) * plot_h
            x = start_x + m_index * (bar_w + bar_gap)
            y = top + plot_h - bar_h
            color = METHOD_COLORS.get(method, "#6b7280")
            parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}" fill="{color}"/>')
            if bar_h > 20:
                parts.append(svg_text(x + bar_w / 2, y - 6, f"{value:.2f}", size=10))

    legend_x = left + 12
    legend_y = height - 34
    for method in methods:
        color = METHOD_COLORS.get(method, "#6b7280")
        parts.append(f'<rect x="{legend_x}" y="{legend_y - 10}" width="12" height="12" fill="{color}"/>')
        parts.append(svg_text(legend_x + 18, legend_y, METHOD_LABELS.get(method, method), size=12, anchor="start"))
        legend_x += 160

    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def genome_length_plot(path: Path, rows: list[dict[str, str]]) -> None:
    width = 900
    height = 500
    left = 88
    right = 28
    top = 62
    bottom = 88
    plot_w = width - left - right
    plot_h = height - top - bottom
    values = {row["condition"]: row for row in rows}
    max_len = max(float(row["max_length"]) for row in rows) if rows else 1.0

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf7"/>',
        svg_text(width / 2, 34, "Genome Lengths By Mutation Condition", size=20, weight="700"),
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#2d3336" stroke-width="1.2"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#2d3336" stroke-width="1.2"/>',
    ]
    for i in range(6):
        value = max_len * i / 5
        y = top + plot_h - (value / max_len) * plot_h
        parts.append(f'<line x1="{left - 5}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" stroke="#d9d4ca" stroke-width="1"/>')
        parts.append(svg_text(left - 10, y + 4, f"{value:.0f}", size=11, anchor="end"))

    group_w = plot_w / len(CONDITIONS)
    for i, condition in enumerate(CONDITIONS):
        row = values.get(condition)
        if not row:
            continue
        x = left + i * group_w + group_w / 2
        min_len = float(row["min_length"])
        mean_len = float(row["mean_length"])
        max_length = float(row["max_length"])
        y_min = top + plot_h - (min_len / max_len) * plot_h
        y_mean = top + plot_h - (mean_len / max_len) * plot_h
        y_max = top + plot_h - (max_length / max_len) * plot_h
        parts.append(f'<line x1="{x:.1f}" y1="{y_max:.1f}" x2="{x:.1f}" y2="{y_min:.1f}" stroke="#2f6f73" stroke-width="3"/>')
        parts.append(f'<circle cx="{x:.1f}" cy="{y_mean:.1f}" r="6" fill="#b24f2a"/>')
        parts.append(svg_text(x, top + plot_h + 32, CONDITION_LABELS[condition], size=12))
        parts.append(svg_text(x, y_mean - 12, f"{mean_len:.0f}", size=10))
    parts.append(svg_text(left + plot_w - 210, height - 28, "Line = min/max, dot = mean", size=12, anchor="start"))
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_summary(path: Path, score_rows: list[dict[str, str]]) -> None:
    lines = ["method\tcondition\trf_distance\tnormalized_rf\tmatching_splits"]
    for row in score_rows:
        lines.append(
            "\t".join(
                [
                    row["method"],
                    row["condition"],
                    row["rf_distance"],
                    row["normalized_rf"],
                    row["matching_splits"],
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    root = project_root()
    plot_dir = ensure_dir(root / "results" / "plots")
    scores = read_tsv(root / "results" / "tree_scores" / "tree_scores.tsv")
    patristic = [row for row in read_tsv(root / "results" / "tree_scores" / "patristic_correlations.tsv") if row.get("patristic_correlation")]
    genome_qc = read_tsv(root / "results" / "simulation_qc" / "phylo_input_qc.tsv")

    if not scores:
        raise FileNotFoundError("No tree score rows found. Run scripts/12_score_trees.py first.")

    bar_chart(
        path=plot_dir / "rf_by_condition_method.svg",
        title="RF Distance By Mutation Condition And Method",
        rows=scores,
        value_key="rf_distance",
        y_label="RF distance",
        y_max=26,
    )
    bar_chart(
        path=plot_dir / "normalized_rf_by_condition_method.svg",
        title="Normalized RF Distance By Mutation Condition And Method",
        rows=scores,
        value_key="normalized_rf",
        y_label="Normalized RF",
        y_max=1.0,
    )
    if patristic:
        bar_chart(
            path=plot_dir / "patristic_correlation_by_condition_method.svg",
            title="Patristic Distance Correlation",
            rows=patristic,
            value_key="patristic_correlation",
            y_label="Correlation",
            y_max=1.0,
        )
    if genome_qc:
        genome_length_plot(plot_dir / "genome_lengths_by_condition.svg", genome_qc)
    write_summary(plot_dir / "score_plot_data.tsv", scores)

    print(f"Wrote plots to {plot_dir}")


if __name__ == "__main__":
    main()
