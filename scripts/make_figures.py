"""
Figure generation for aevol_phylo_mutrate pilot benchmark report.
Produces four figures from existing results data.

Usage:
    python scripts/make_figures.py

Output: results/plots/fig1_genome_lengths.png
                       fig2_rf_distances.png
                       fig3_true_tree.png
                       fig4_patristic_correlations.png
"""

import os
import sys

# Always resolve paths relative to the aevol_phylo_mutrate project root,
# regardless of where the script is invoked from.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
os.chdir(PROJECT_ROOT)

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

OUT_DIR = "results/plots"
os.makedirs(OUT_DIR, exist_ok=True)

# ── Shared style ──────────────────────────────────────────────────────────────
CONDITIONS = ["very_low", "low", "normal", "high", "very_high"]
CONDITION_LABELS = ["very_low\n(0.0625x)", "low\n(0.25x)", "normal\n(1x)", "high\n(4x)", "very_high\n(16x)"]
METHOD_COLORS = {
    "mash_bionj":   "#E07B39",
    "mafft_iqtree": "#3A7EBF",
    "mauve_iqtree":  "#4CAF70",
}
METHOD_LABELS = {
    "mash_bionj":   "Mash/BioNJ",
    "mafft_iqtree": "MAFFT/IQ-TREE",
    "mauve_iqtree":  "Mauve/IQ-TREE",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
})

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1: Genome length distribution across mutation-rate conditions
# ─────────────────────────────────────────────────────────────────────────────
lengths = pd.read_csv("results/simulation_qc/genome_lengths.tsv", sep="\t")
lengths["condition"] = pd.Categorical(lengths["condition"], categories=CONDITIONS, ordered=True)

fig, ax = plt.subplots(figsize=(9, 5.5))

data_by_cond = [lengths[lengths["condition"] == c]["length"].values for c in CONDITIONS]
bp = ax.boxplot(
    data_by_cond,
    positions=range(len(CONDITIONS)),
    patch_artist=True,
    widths=0.5,
    medianprops={"color": "black", "linewidth": 2},
    whiskerprops={"linewidth": 1.2},
    capprops={"linewidth": 1.2},
    flierprops={"marker": "o", "markersize": 5, "alpha": 0.7},
)
box_colors = ["#D0E8FF", "#A8D0FF", "#7CB9FF", "#FF9966", "#FF5533"]
for patch, color in zip(bp["boxes"], box_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.85)

# Overlay individual points with jitter
for i, cond in enumerate(CONDITIONS):
    vals = lengths[lengths["condition"] == cond]["length"].values
    jitter = np.random.default_rng(seed=42).uniform(-0.15, 0.15, size=len(vals))
    ax.scatter(i + jitter, vals, color="black", alpha=0.5, s=22, zorder=3)

ax.set_xticks(range(len(CONDITIONS)))
ax.set_xticklabels(CONDITION_LABELS, fontsize=10)
ax.set_ylabel("Genome length (bp)", fontsize=11)
ax.set_xlabel("Mutation-rate condition", fontsize=11)
ax.set_title(
    "Figure 1: Final genome length by mutation-rate condition (rep_001)",
    fontweight="bold", pad=12, fontsize=12,
)

# Annotate means inside each box (at median height), avoiding overlap with axis labels
for i, (cond, box) in enumerate(zip(CONDITIONS, bp["boxes"])):
    mean = lengths[lengths["condition"] == cond]["length"].mean()
    # Place label just above the median line
    median_y = lengths[lengths["condition"] == cond]["length"].median()
    ax.text(
        i, median_y + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.025,
        f"mean={mean:.0f}",
        ha="center", va="bottom", fontsize=8, color="#222",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.6, pad=1),
    )

plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig1_genome_lengths.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved fig1_genome_lengths.png")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: Normalized RF distance by method and condition
# ─────────────────────────────────────────────────────────────────────────────
scores = pd.read_csv("results/tree_scores/tree_scores.tsv", sep="\t")
scores["condition"] = pd.Categorical(scores["condition"], categories=CONDITIONS, ordered=True)

fig, ax = plt.subplots(figsize=(9, 5))

methods = list(METHOD_COLORS.keys())
n_methods = len(methods)
width = 0.25
x = np.arange(len(CONDITIONS))

for i, method in enumerate(methods):
    sub = scores[scores["method"] == method].sort_values("condition")
    sub = sub.set_index("condition").reindex(CONDITIONS)
    vals = sub["normalized_rf"].values
    offset = (i - 1) * width
    bars = ax.bar(
        x + offset, vals,
        width=width,
        color=METHOD_COLORS[method],
        label=METHOD_LABELS[method],
        alpha=0.85,
        edgecolor="white",
        linewidth=0.5,
    )
    # Mark missing bars
    for j, (v, cond) in enumerate(zip(vals, CONDITIONS)):
        if pd.isna(v):
            ax.text(
                x[j] + offset, 0.02, "n/a",
                ha="center", va="bottom", fontsize=8,
                color=METHOD_COLORS[method], fontweight="bold",
            )

ax.set_xticks(x)
ax.set_xticklabels(CONDITION_LABELS)
ax.set_ylabel("Normalized Robinson-Foulds distance")
ax.set_xlabel("Mutation-rate condition")
ax.set_ylim(0, 1.0)
ax.set_title("Figure 2: Tree recovery (normalized RF distance) by method and condition\n(rep_001; lower = better; 0 = perfect recovery)", fontweight="bold", pad=10)
ax.axhline(0, color="black", linewidth=0.6, linestyle="-")
ax.legend(loc="upper left", framealpha=0.9)

# Add value labels on top of bars
for i, method in enumerate(methods):
    sub = scores[scores["method"] == method].sort_values("condition")
    sub = sub.set_index("condition").reindex(CONDITIONS)
    vals = sub["normalized_rf"].values
    offset = (i - 1) * width
    for j, v in enumerate(vals):
        if not pd.isna(v) and v > 0:
            ax.text(x[j] + offset, v + 0.01, f"{v:.3f}", ha="center", va="bottom", fontsize=7.5)

plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig2_rf_distances.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved fig2_rf_distances.png")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 3: True 16-tip balanced tree drawn as a cladogram
# ─────────────────────────────────────────────────────────────────────────────
# Tree topology (from true_tree_16tips.nwk):
#   ((((S01,S02),(S03,S04)),((S05,S06),(S07,S08))),
#    (((S09,S10),(S11,S12)),((S13,S14),(S15,S16))))
#
# Layout: tips get y = 0..15 (top-to-bottom), x = their depth level (0=root, 4=tips).
# Each internal node's y = mean of its children's y values.
# Lines: horizontal from parent to node, then vertical between children.

def clade_layout(labels, depth=0):
    """
    Recursively compute (node_x, node_y) for every node and collect line segments.
    Returns: (node_y, tip_positions, segments)
      tip_positions: {label: y}  (filled as a side-effect via closure)
      segments: list of ((x0,y0),(x1,y1))
    """
    tip_positions = {}
    segments = []
    _tip_counter = [0]  # mutable counter via list

    def recurse(label_list, depth):
        if len(label_list) == 1:
            y = _tip_counter[0]
            _tip_counter[0] += 1
            tip_positions[label_list[0]] = y
            return depth, y           # (node_x, node_y)

        mid = len(label_list) // 2
        left_labels  = label_list[:mid]
        right_labels = label_list[mid:]

        lx, ly = recurse(left_labels,  depth + 1)
        rx, ry = recurse(right_labels, depth + 1)

        node_x = depth
        node_y = (ly + ry) / 2

        # Horizontal lines: parent node → each child node
        segments.append(((node_x, node_y), (lx, ly)))
        segments.append(((node_x, node_y), (rx, ry)))
        # Vertical connector at this node's x
        segments.append(((node_x, ly), (node_x, ry)))

        return node_x, node_y

    recurse(labels, 0)
    return tip_positions, segments

labels_ordered = [f"S{i:02d}" for i in range(1, 17)]
tip_positions, segments = clade_layout(labels_ordered)

fig, ax = plt.subplots(figsize=(7, 7))
for (x0, y0), (x1, y1) in segments:
    ax.plot([x0, x1], [y0, y1], color="#2255AA", linewidth=1.8, solid_capstyle="round")

for tip, y in tip_positions.items():
    ax.text(4.08, y, tip, va="center", ha="left", fontsize=9)

ax.set_xlim(-0.3, 5.2)
ax.set_ylim(-1, 16)
ax.axis("off")
ax.set_title(
    "Figure 3: Known true 16-tip balanced species tree\n"
    "(all branch lengths equal; used as ground truth in all conditions)",
    fontweight="bold", pad=10,
)

for d, label in enumerate(["Root", "Depth 1", "Depth 2", "Depth 3", "Tips"]):
    ax.text(d, -0.8, label, ha="center", fontsize=8, color="#666", style="italic")

plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig3_true_tree.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved fig3_true_tree.png")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4: Patristic correlation (secondary diagnostic) + RF overlay line plot
# ─────────────────────────────────────────────────────────────────────────────
patristic = pd.read_csv("results/tree_scores/patristic_correlations.tsv", sep="\t")
patristic["patristic_correlation"] = pd.to_numeric(patristic["patristic_correlation"], errors="coerce")
patristic["condition"] = pd.Categorical(patristic["condition"], categories=CONDITIONS, ordered=True)

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

# Left panel: RF distance (line plot for trends)
ax = axes[0]
for method in methods:
    sub = scores[scores["method"] == method].sort_values("condition")
    sub = sub.set_index("condition").reindex(CONDITIONS)
    x_vals = range(len(CONDITIONS))
    y_vals = sub["normalized_rf"].values
    ax.plot(x_vals, y_vals, marker="o", color=METHOD_COLORS[method],
            label=METHOD_LABELS[method], linewidth=2, markersize=7)
    # Mark missing
    for j, v in enumerate(y_vals):
        if pd.isna(v):
            ax.scatter(j, 0.5, marker="x", s=80, color=METHOD_COLORS[method], zorder=5)
            ax.text(j, 0.52, "n/a", ha="center", fontsize=8, color=METHOD_COLORS[method])

ax.set_xticks(range(len(CONDITIONS)))
ax.set_xticklabels(CONDITION_LABELS, fontsize=9)
ax.set_ylabel("Normalized RF distance (lower = better)")
ax.set_ylim(-0.05, 1.0)
ax.set_title("(a) RF distance trends", fontweight="bold")
ax.legend(fontsize=9)

# Right panel: patristic correlation
ax = axes[1]
for method in methods:
    sub = patristic[patristic["method"] == method].sort_values("condition")
    sub = sub.set_index("condition").reindex(CONDITIONS)
    x_vals = range(len(CONDITIONS))
    y_vals = sub["patristic_correlation"].values
    ax.plot(x_vals, y_vals, marker="s", color=METHOD_COLORS[method],
            label=METHOD_LABELS[method], linewidth=2, markersize=7, linestyle="--")
    for j, v in enumerate(y_vals):
        if pd.isna(v):
            ax.scatter(j, 0.5, marker="x", s=80, color=METHOD_COLORS[method], zorder=5)
            ax.text(j, 0.52, "n/a", ha="center", fontsize=8, color=METHOD_COLORS[method])

ax.set_xticks(range(len(CONDITIONS)))
ax.set_xticklabels(CONDITION_LABELS, fontsize=9)
ax.set_ylabel("Patristic distance correlation (higher = better)")
ax.set_ylim(0, 1.05)
ax.set_title("(b) Patristic distance correlations (secondary)", fontweight="bold")
ax.legend(fontsize=9)

fig.suptitle("Figure 4: RF distance trends and patristic distance correlations by method (rep_001)",
             fontweight="bold", y=1.02)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig4_summary_trends.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved fig4_summary_trends.png")

print("\nAll figures written to", OUT_DIR)
