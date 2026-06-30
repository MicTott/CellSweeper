# =============================================================================
# Benchmark: Where does clusterCellTypes spend its time?
# =============================================================================
#
# Runs each step of clusterCellTypes individually on the spatialDLPFC
# snRNA-seq dataset (~77k cells) and reports wall-clock time per step.
#
# Usage: source this script in an R session with CellSweeper loaded.
# =============================================================================

devtools::load_all()
library(spatialLIBD)
library(SpatialExperiment)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle)
library(BiocNeighbors)
library(BiocSingular)

# =============================================================================
# 1. Load data
# =============================================================================

cat("=== Loading data ===\n")
t0 <- system.time({
    zip_path <- fetch_data("spatialDLPFC_snRNAseq")
    tmpdir <- tempdir()
    unzip(zip_path, exdir = tmpdir)
    sce <- HDF5Array::loadHDF5SummarizedExperiment(
        file.path(tmpdir, "sce_DLPFC_annotated")
    )
})
cat("  Load + unzip:", t0["elapsed"], "s\n")

cat("  Dimensions:", dim(sce), "\n")
cat("  Cell types:\n")
print(table(sce$cellType_broad_hc))

# Convert to SPE
cat("\n=== Converting to SpatialExperiment ===\n")
t1 <- system.time({
    counts_mat <- as(counts(sce), "dgCMatrix")
    coords <- matrix(rnorm(ncol(sce) * 2), ncol = 2,
                     dimnames = list(colnames(sce), c("x", "y")))
    spe <- SpatialExperiment(
        assays = list(counts = counts_mat),
        spatialCoords = coords,
        colData = colData(sce)
    )
    spe <- addPerCellQCMetrics(spe,
        subsets = list(mito = grep("^MT-", rownames(spe))))
})
cat("  Realize + build SPE + QC:", t1["elapsed"], "s\n")
cat("  Final dims:", dim(spe), "\n\n")

# Store expert labels
spe$expert_broad <- spe$cellType_broad_hc

# =============================================================================
# 2. Benchmark each step of clusterCellTypes (no regression)
# =============================================================================

cat("====================================================\n")
cat("STEP-BY-STEP BENCHMARK (no covariate regression)\n")
cat("====================================================\n\n")

# --- Gene filtering ---
cat("--- Gene filtering ---\n")
t_filter <- system.time({
    genes_to_use <- rownames(spe)
    mito_genes <- grep("^MT-|^mt-", genes_to_use, value = TRUE)
    genes_to_use <- setdiff(genes_to_use, mito_genes)
    ribo_genes <- grep("^RP[SL]|^Rp[sl]", genes_to_use, value = TRUE)
    genes_to_use <- setdiff(genes_to_use, ribo_genes)
    spe_sub <- spe[genes_to_use, ]
})
cat("  Excluded", length(mito_genes), "mito +", length(ribo_genes), "ribo genes\n")
cat("  Remaining:", nrow(spe_sub), "genes\n")
cat("  Time:", t_filter["elapsed"], "s\n\n")

# --- Normalization ---
cat("--- Log-normalization ---\n")
t_norm <- system.time({
    spe_sub <- scater::logNormCounts(spe_sub)
})
cat("  Time:", t_norm["elapsed"], "s\n\n")

# --- Feature selection (HVGs) ---
cat("--- Feature selection (modelGeneVar + top 2000 HVGs) ---\n")
t_hvg <- system.time({
    dec <- scran::modelGeneVar(spe_sub, assay.type = "logcounts")
    hvgs <- scran::getTopHVGs(dec, n = 2000)
})
cat("  Selected", length(hvgs), "HVGs\n")
cat("  Time:", t_hvg["elapsed"], "s\n\n")

# --- PCA (exact) ---
cat("--- PCA (exact SVD, 30 PCs) ---\n")
t_pca_exact <- system.time({
    spe_sub_exact <- scater::runPCA(spe_sub, subset_row = hvgs,
                                     ncomponents = 30,
                                     exprs_values = "logcounts")
})
cat("  Time:", t_pca_exact["elapsed"], "s\n\n")

# --- PCA (IRLBA) ---
cat("--- PCA (IrlbaParam, 30 PCs) ---\n")
t_pca_irlba <- system.time({
    spe_sub_irlba <- scater::runPCA(spe_sub, subset_row = hvgs,
                                     ncomponents = 30,
                                     exprs_values = "logcounts",
                                     BSPARAM = IrlbaParam())
})
cat("  Time:", t_pca_irlba["elapsed"], "s\n\n")

# Use IRLBA result going forward
spe_sub <- spe_sub_irlba

# --- Leverage score computation ---
cat("--- Leverage score computation ---\n")
t_leverage <- system.time({
    pca_mat <- reducedDim(spe_sub, "PCA")
    col_sd <- apply(pca_mat, 2, sd)
    col_sd[col_sd == 0] <- 1
    pca_scaled <- sweep(pca_mat, 2, col_sd, "/")
    leverage <- rowSums(pca_scaled^2)
    lev_prob <- leverage / sum(leverage)
})
cat("  Time:", t_leverage["elapsed"], "s\n\n")

# --- Subsampling ---
cat("--- Subsampling (20k cells) ---\n")
set.seed(42)
t_subsample <- system.time({
    sub_idx <- sample(ncol(spe_sub), 20000, prob = lev_prob)
    spe_sketch <- spe_sub[, sub_idx]
})
cat("  Time:", t_subsample["elapsed"], "s\n\n")

# --- Clustering: Leiden on subsample ---
cat("--- Leiden clustering on 20k subsample (k=25, res=0.5) ---\n")
t_leiden_sub <- system.time({
    blusparam_sub <- bluster::NNGraphParam(
        cluster.fun = "leiden",
        cluster.args = list(resolution_parameter = 0.5),
        k = 25
    )
    clusters_sub <- scran::clusterCells(spe_sketch,
        use.dimred = "PCA", BLUSPARAM = blusparam_sub)
})
cat("  Found", length(unique(clusters_sub)), "clusters\n")
cat("  Time:", t_leiden_sub["elapsed"], "s\n\n")

# --- Leiden clustering on FULL dataset ---
cat("--- Leiden clustering on FULL dataset (", ncol(spe_sub), "cells) ---\n")
t_leiden_full <- system.time({
    blusparam_full <- bluster::NNGraphParam(
        cluster.fun = "leiden",
        cluster.args = list(resolution_parameter = 0.5),
        k = 25
    )
    clusters_full <- scran::clusterCells(spe_sub,
        use.dimred = "PCA", BLUSPARAM = blusparam_full)
})
cat("  Found", length(unique(clusters_full)), "clusters\n")
cat("  Time:", t_leiden_full["elapsed"], "s\n\n")

# --- Leiden on full with HNSW ---
cat("--- Leiden on FULL dataset with HnswParam ---\n")
t_leiden_hnsw <- system.time({
    blusparam_hnsw <- bluster::NNGraphParam(
        cluster.fun = "leiden",
        cluster.args = list(resolution_parameter = 0.5),
        k = 25,
        BNPARAM = HnswParam()
    )
    clusters_hnsw <- scran::clusterCells(spe_sub,
        use.dimred = "PCA", BLUSPARAM = blusparam_hnsw)
})
cat("  Found", length(unique(clusters_hnsw)), "clusters\n")
cat("  Time:", t_leiden_hnsw["elapsed"], "s\n\n")

# --- Label projection (centroid assignment) ---
cat("--- Label projection (centroid nearest-neighbor for remaining cells) ---\n")
t_project <- system.time({
    pca_all <- reducedDim(spe_sub, "PCA")
    pca_sk <- pca_all[sub_idx, , drop = FALSE]
    cls_labels <- as.character(clusters_sub)
    unique_cls <- unique(cls_labels)
    centroids <- do.call(rbind, lapply(unique_cls, function(cl) {
        colMeans(pca_sk[cls_labels == cl, , drop = FALSE])
    }))
    nn_result <- BiocNeighbors::queryKNN(
        X = centroids, query = pca_all, k = 1)
    projected_labels <- unique_cls[nn_result$index[, 1]]
})
cat("  Time:", t_project["elapsed"], "s\n\n")

# =============================================================================
# 3. Summary
# =============================================================================

cat("====================================================\n")
cat("SUMMARY (seconds)\n")
cat("====================================================\n")
timings <- data.frame(
    Step = c(
        "Gene filtering",
        "Log-normalization",
        "Feature selection (HVGs)",
        "PCA (exact SVD)",
        "PCA (IRLBA)",
        "Leverage scores",
        "Subsampling (20k)",
        "Leiden on 20k subsample",
        "Leiden on FULL dataset",
        "Leiden on FULL + HNSW",
        "Label projection"
    ),
    Seconds = c(
        t_filter["elapsed"],
        t_norm["elapsed"],
        t_hvg["elapsed"],
        t_pca_exact["elapsed"],
        t_pca_irlba["elapsed"],
        t_leverage["elapsed"],
        t_subsample["elapsed"],
        t_leiden_sub["elapsed"],
        t_leiden_full["elapsed"],
        t_leiden_hnsw["elapsed"],
        t_project["elapsed"]
    )
)
print(timings, row.names = FALSE)

cat("\n--- Sketch pipeline total (IRLBA + leverage + subsample + Leiden 20k + project) ---\n")
sketch_total <- t_filter["elapsed"] + t_norm["elapsed"] + t_hvg["elapsed"] +
    t_pca_irlba["elapsed"] + t_leverage["elapsed"] + t_subsample["elapsed"] +
    t_leiden_sub["elapsed"] + t_project["elapsed"]
cat("  ", round(sketch_total, 1), "s\n")

cat("\n--- Full pipeline total (exact PCA + Leiden full) ---\n")
full_total <- t_filter["elapsed"] + t_norm["elapsed"] + t_hvg["elapsed"] +
    t_pca_exact["elapsed"] + t_leiden_full["elapsed"]
cat("  ", round(full_total, 1), "s\n")

cat("\n--- Speedup (full / sketch): ", round(full_total / sketch_total, 1), "x ---\n")

# =============================================================================
# 4. Compare cluster quality vs expert annotations
# =============================================================================

cat("\n====================================================\n")
cat("CLUSTER QUALITY vs EXPERT ANNOTATIONS\n")
cat("====================================================\n\n")

expert <- spe$expert_broad

cat("--- Full Leiden (all cells) ---\n")
print(table(clusters_full, expert))

cat("\n--- Sketch Leiden (leverage 20k, projected) ---\n")
print(table(projected_labels, expert))

cat("\nBenchmark complete.\n")
