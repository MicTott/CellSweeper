# =============================================================================
# CellSweeper Demo Script
# =============================================================================
#
# This script walks through the full CellSweeper pipeline using synthetic data
# with 3 spatially separated cell-type clusters and injected outliers.
#
# It demonstrates:
#   1. Full pipeline (runCellSweeper)
#   2. Step-by-step workflow
#   3. Performance options (subsample, leverage, approximate NN, etc.)
#   4. Visualization
#
# No external data download required - everything runs on synthetic data.
# =============================================================================

# devtools::load_all() for development; use library(CellSweeper) if installed
devtools::load_all()
library(SpatialExperiment)
library(SummarizedExperiment)
library(scuttle)

# =============================================================================
# 1. Create synthetic data
# =============================================================================
#
# Build a SpatialExperiment with:
#   - 3 spatially separated clusters (blobs in 2D tissue space)
#   - 200 genes (including 10 MT- genes, cluster-specific markers)
#   - Injected outliers: near-zero library size + morphology anomalies
#   - QC metrics already computed (sum, detected, subsets_mito_percent)
#   - Morphology columns (Area_um, log2AspectRatio, log2CountArea)

set.seed(42)
n_cells_per_cluster <- 200
n_genes <- 200
n_outliers <- 10
n_clusters <- 3
n_cells <- n_cells_per_cluster * n_clusters

# Spatial coordinates: 3 well-separated blobs
centers <- matrix(c(0, 0, 100, 0, 50, 100), ncol = 2, byrow = TRUE)
coords <- do.call(rbind, lapply(seq_len(n_clusters), function(i) {
    matrix(rnorm(n_cells_per_cluster * 2, mean = 0, sd = 10), ncol = 2) +
        matrix(rep(centers[i, ], each = n_cells_per_cluster), ncol = 2)
}))
colnames(coords) <- c("x_centroid", "y_centroid")

# Gene expression with cluster-specific markers
cluster_labels <- rep(paste0("cluster_", 1:3), each = n_cells_per_cluster)
base_means <- c(5, 3, 1)
counts_matrix <- do.call(cbind, lapply(seq_len(n_clusters), function(i) {
    matrix(rpois(n_genes * n_cells_per_cluster, lambda = base_means[i]),
           nrow = n_genes)
}))
for (i in 1:n_clusters) {
    gene_idx <- ((i - 1) * 10 + 1):(i * 10)
    col_idx <- ((i - 1) * n_cells_per_cluster + 1):(i * n_cells_per_cluster)
    counts_matrix[gene_idx, col_idx] <- counts_matrix[gene_idx, col_idx] + 10
}

# Name genes (last 10 are mitochondrial)
gene_names <- paste0("Gene", 1:n_genes)
gene_names[(n_genes - 9):n_genes] <- paste0("MT-", gene_names[(n_genes - 9):n_genes])
rownames(counts_matrix) <- gene_names
colnames(counts_matrix) <- paste0("cell_", 1:n_cells)

# Inject outliers: near-zero library size
outlier_flags <- rep(FALSE, n_cells)
for (i in 1:n_clusters) {
    start <- (i - 1) * n_cells_per_cluster + 1
    idx <- start:(start + n_outliers - 1)
    outlier_flags[idx] <- TRUE
    counts_matrix[, idx] <- 0
    counts_matrix[1:2, idx] <- 1
}

# Build SpatialExperiment
spe <- SpatialExperiment(
    assays = list(counts = as(counts_matrix, "dgCMatrix")),
    spatialCoords = coords,
    sample_id = rep("sample1", n_cells)
)
spe <- addPerCellQCMetrics(spe,
    subsets = list(mito = grep("^MT-", rownames(spe))))

# Add morphology columns
colData(spe)$Area_um <- pmax(50 + spe$sum * 0.05 + rnorm(n_cells, 0, 10), 10)
colData(spe)$log2AspectRatio <- rnorm(n_cells, 0, 0.5)
colData(spe)$CountArea <- spe$sum / colData(spe)$Area_um
colData(spe)$log2CountArea <- log2(colData(spe)$CountArea + 1)
for (i in 1:n_clusters) {
    start <- (i - 1) * n_cells_per_cluster + 1
    idx <- start:(start + n_outliers - 1)
    colData(spe)$Area_um[idx] <- 5
    colData(spe)$log2AspectRatio[idx] <- 3.0
}

# Store ground truth
spe$true_outlier <- outlier_flags
spe$true_cluster <- cluster_labels
cat("Synthetic data:", ncol(spe), "cells,", nrow(spe), "genes\n")
cat("True clusters:", paste(unique(spe$true_cluster), collapse = ", "), "\n")
cat("Injected outliers:", sum(spe$true_outlier), "\n\n")

# Peek at colData
head(colData(spe)[, c("sum", "detected", "subsets_mito_percent",
                       "Area_um", "true_cluster", "true_outlier")])


# =============================================================================
# 2. Full pipeline with runCellSweeper
# =============================================================================
#
# The easiest way to run CellSweeper: one function call.
# Here we use pre-existing cluster labels for speed.

spe_full <- runCellSweeper(spe,
    use_existing = "true_cluster",    # skip clustering, use known labels
    morpho_metrics = c("Area_um", "log2AspectRatio", "log2CountArea"),
    k = 30,
    z_threshold = 3,
    n_permutations = 100,
    verbose = TRUE
)

cat("\n--- Results ---\n")
cat("Artifact clusters:", sum(spe_full$cluster_artifact_flag), "cells\n")
cat("Spatial outliers: ", sum(spe_full$cellsweeper_outlier), "cells\n")
cat("Total flagged:    ", sum(spe_full$cellsweeper_qc), "cells\n")

# Compare to ground truth
cat("\nGround truth outliers recovered:\n")
print(table(
    flagged = spe_full$cellsweeper_outlier,
    true_outlier = spe_full$true_outlier
))


# =============================================================================
# 3. Step-by-step workflow
# =============================================================================

# --- Level 1: Global pre-filtering ---
spe2 <- globalFilter(spe, min_counts = 5, min_genes = 5)
cat("\nAfter global filter:", ncol(spe2), "cells remain\n")

# --- Level 2a: Clustering ---
# Default: Leiden on SNN graph, with QC covariate regression
spe2 <- clusterCellTypes(spe2,
    resolution = 0.5,
    regress_covariates = FALSE,  # faster for demo
    exclude_mito = TRUE,
    exclude_ribo = TRUE
)
cat("Clusters found:", length(unique(spe2$cell_cluster)), "\n")
print(table(spe2$cell_cluster))

# --- Level 2b+c: Flag artifact clusters ---
spe2 <- flagArtifactClusters(spe2,
    metrics = c("sum", "detected"),
    n_permutations = 100,
    seed = 42
)
cat("\nArtifact clusters flagged:", sum(spe2$cluster_artifact_flag), "cells\n")
cat("Per-cluster spatial homogeneity:\n")
print(data.frame(
    cluster = unique(spe2$cell_cluster),
    homogeneity = tapply(spe2$cluster_spatial_homogeneity,
                         spe2$cell_cluster, unique)
))

# --- Level 3: Within-cluster spatial QC ---
spe2 <- clusterLocalOutliers(spe2,
    metrics = c("sum", "detected"),
    morpho_metrics = c("Area_um", "log2AspectRatio", "log2CountArea"),
    k = 30,
    z_threshold = 3,
    combine = "union"
)
cat("\nSpatial outliers:", sum(spe2$cellsweeper_outlier), "cells\n")

# Final combined flag
spe2$cellsweeper_qc <- spe2$cluster_artifact_flag | spe2$cellsweeper_outlier
spe2$cellsweeper_qc[is.na(spe2$cellsweeper_qc)] <- FALSE
cat("Total flagged:", sum(spe2$cellsweeper_qc), "/", ncol(spe2), "\n")


# =============================================================================
# 4. Performance options
# =============================================================================
#
# For large datasets (>100k cells), several speedups are available.
# All are opt-in; defaults produce identical results to the base pipeline.

# --- 4a. Sketch-based clustering with leverage score sampling ---
# Subsample 200 cells, cluster them, project labels to the rest.
# Leverage sampling preferentially selects rare/outlying cells.
spe_sketch <- clusterCellTypes(spe,
    subsample = 200,                    # sketch: cluster 200, project to all
    subsample_method = "leverage",      # default: leverage score sampling
    regress_covariates = FALSE
)
cat("\nSketch clustering (leverage):",
    length(unique(spe_sketch$cell_cluster)), "cluster(s)\n")

# With alpha mixing (50% leverage + 50% uniform)
spe_mixed <- clusterCellTypes(spe,
    subsample = 200,
    leverage_alpha = 0.5,
    regress_covariates = FALSE
)

# Pure random sampling for comparison
spe_random <- clusterCellTypes(spe,
    subsample = 200,
    subsample_method = "random",
    regress_covariates = FALSE
)

# --- 4b. Approximate nearest neighbors ---
# Use HNSW (Hierarchical Navigable Small World) for faster kNN
# in clustering, outlier detection, and spatial homogeneity.
spe_approx <- clusterCellTypes(spe,
    BNPARAM = BiocNeighbors::HnswParam(),
    regress_covariates = FALSE
)
cat("Approximate NN clustering:",
    length(unique(spe_approx$cell_cluster)), "cluster(s)\n")

# --- 4c. Combined: sketch + approximate NN + fast PCA ---
# This is the "maximum speed" configuration for very large datasets.
if (requireNamespace("BiocSingular", quietly = TRUE)) {
    spe_fast <- clusterCellTypes(spe,
        subsample = 200,
        subsample_method = "leverage",
        BNPARAM = BiocNeighbors::HnswParam(),
        BSPARAM = BiocSingular::IrlbaParam(),
        regress_covariates = FALSE
    )
    cat("Fast config:",
        length(unique(spe_fast$cell_cluster)), "cluster(s)\n")
}


# =============================================================================
# 5. Visualization
# =============================================================================

if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)

    # Use the step-by-step result (spe2) for visualization
    # since it has all QC columns populated

    # --- Spatial outlier plots ---
    p1 <- plotSpatialOutliers(spe2, metric = "sum")
    print(p1)

    p2 <- plotSpatialOutliers(spe2, metric = "detected")
    print(p2)

    # --- Cluster QC distributions ---
    p3 <- plotClusterQC(spe2, metric = "sum")
    print(p3)

    p4 <- plotClusterQC(spe2, metric = "Area_um")
    print(p4)

    # --- Spatial dispersion ---
    p5 <- plotSpatialDispersion(spe2)
    print(p5)

    # --- Cluster summary ---
    summary_plots <- plotClusterSummary(spe2)
    print(summary_plots$mahalanobis)
    print(summary_plots$homogeneity)

    cat("\nAll plots generated successfully.\n")
}


# =============================================================================
# 6. Inspect per-metric results
# =============================================================================

# Z-scores and outlier flags are stored per-metric in colData
zscore_cols <- grep("_zscore$", colnames(colData(spe2)), value = TRUE)
outlier_cols <- grep("_outlier$", colnames(colData(spe2)), value = TRUE)

cat("\nZ-score columns:", paste(zscore_cols, collapse = ", "), "\n")
cat("Outlier columns:", paste(outlier_cols, collapse = ", "), "\n")

# Distribution of z-scores
for (col in zscore_cols) {
    vals <- colData(spe2)[[col]]
    cat(sprintf("  %s: mean=%.2f, sd=%.2f, range=[%.1f, %.1f]\n",
                col, mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE),
                min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))
}

# Per-metric outlier counts
for (col in outlier_cols) {
    n <- sum(colData(spe2)[[col]], na.rm = TRUE)
    cat(sprintf("  %s: %d flagged\n", col, n))
}

cat("\nDemo complete.\n")
