from __future__ import annotations

import argparse
from dataclasses import dataclass

from pipeline_common import ensure_dir, project_root, read_run_config, rel_path, write_tsv


@dataclass
class Edge:
    edge_id: str
    parent: str
    child: str
    depth: int
    branch_generations: int
    is_tip: bool


class BalancedTreeBuilder:
    def __init__(self, tips: int, branch_generations: int):
        if tips < 2 or tips & (tips - 1):
            raise ValueError("This pilot scaffold expects tips to be a power of two.")
        self.tips = [f"S{i:02d}" for i in range(1, tips + 1)]
        self.branch_generations = branch_generations
        self.internal_counter = 0
        self.edge_counter = 0
        self.edges: list[Edge] = []

    def new_internal(self) -> str:
        self.internal_counter += 1
        return f"N{self.internal_counter:02d}"

    def add_edge(self, parent: str, child: str, depth: int) -> None:
        self.edge_counter += 1
        self.edges.append(
            Edge(
                edge_id=f"E{self.edge_counter:03d}",
                parent=parent,
                child=child,
                depth=depth,
                branch_generations=self.branch_generations,
                is_tip=child.startswith("S"),
            )
        )

    def build_child(self, parent: str, child_tips: list[str], depth: int) -> str:
        if len(child_tips) == 1:
            child = child_tips[0]
            self.add_edge(parent, child, depth)
            return f"{child}:{self.branch_generations}"

        child = self.new_internal()
        self.add_edge(parent, child, depth)
        return self.build_subtree(child, child_tips, depth)

    def build_subtree(self, node: str, node_tips: list[str], depth: int) -> str:
        midpoint = len(node_tips) // 2
        left = self.build_child(node, node_tips[:midpoint], depth + 1)
        right = self.build_child(node, node_tips[midpoint:], depth + 1)
        if node == "root":
            return f"({left},{right}){node};"
        return f"({left},{right}){node}:{self.branch_generations}"

    def build(self) -> tuple[str, list[Edge]]:
        newick = self.build_subtree("root", self.tips, 0)
        return newick, self.edges


def parse_args() -> argparse.Namespace:
    root = project_root()
    config = read_run_config(root)
    parser = argparse.ArgumentParser(description="Create a balanced true species tree.")
    parser.add_argument("--tips", type=int, default=int(config["tips"]))
    parser.add_argument(
        "--branch-generations",
        type=int,
        default=int(config["branch_generations"]),
    )
    parser.add_argument("--out-tree", default=str(config["true_tree"]))
    parser.add_argument("--out-edges", default=str(config["tree_edges"]))
    return parser.parse_args()


def main() -> None:
    root = project_root()
    args = parse_args()
    builder = BalancedTreeBuilder(args.tips, args.branch_generations)
    newick, edges = builder.build()

    tree_path = rel_path(root, args.out_tree)
    edge_path = rel_path(root, args.out_edges)
    ensure_dir(tree_path.parent)
    tree_path.write_text(newick + "\n", encoding="utf-8")

    rows = [
        {
            "edge_id": edge.edge_id,
            "parent": edge.parent,
            "child": edge.child,
            "depth": edge.depth,
            "branch_generations": edge.branch_generations,
            "is_tip": str(edge.is_tip).lower(),
        }
        for edge in edges
    ]
    write_tsv(
        edge_path,
        rows,
        ["edge_id", "parent", "child", "depth", "branch_generations", "is_tip"],
    )

    print(f"Wrote {tree_path}")
    print(f"Wrote {edge_path} with {len(rows)} branches")


if __name__ == "__main__":
    main()

