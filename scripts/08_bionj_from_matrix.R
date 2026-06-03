#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 08_bionj_from_matrix.R distance_matrix.tsv inferred_tree.nwk", call. = FALSE)
}

matrix_path <- args[[1]]
tree_path <- args[[2]]

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("R package 'ape' is required for ape::bionj().", call. = FALSE)
}

dist_table <- utils::read.table(
  matrix_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

dist_matrix <- as.matrix(dist_table)
storage.mode(dist_matrix) <- "double"

if (!identical(rownames(dist_matrix), colnames(dist_matrix))) {
  stop("Distance matrix row and column labels do not match.", call. = FALSE)
}
if (any(!is.finite(dist_matrix))) {
  stop("Distance matrix contains non-finite values.", call. = FALSE)
}
if (max(abs(dist_matrix - t(dist_matrix))) > 1e-10) {
  stop("Distance matrix is not symmetric.", call. = FALSE)
}
diag(dist_matrix) <- 0

tree <- ape::bionj(stats::as.dist(dist_matrix))
dir.create(dirname(tree_path), recursive = TRUE, showWarnings = FALSE)
ape::write.tree(tree, file = tree_path)

cat("Wrote", tree_path, "\n")
