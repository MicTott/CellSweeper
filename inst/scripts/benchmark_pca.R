# =============================================================================
# Benchmark: Why is PCA so slow? What can we do about it?
# =============================================================================
#
# Investigates why IRLBA isn't faster than exact SVD on the DLPFC dataset
# and tests alternative approaches.
# =============================================================================

devtools::load_all()
library(spatialLIBD)
library(SpatialExperiment)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle)
library(BiocSingular)
library(Matrix)

# =============================================================================
# 1. Load and prepare data (reuse from benchmark_clustering.R)
# =============================================================================

cat("=== Loading data ===\n")
zip_path <- fetch_data("spatialDLPFC_snRNAseq")
tmpdir <- tempdir()
unzip(zip_path, exdir = tmpdir)
sce <- HDF5Array::loadHDF5SummarizedExperiment(
    file.path(tmpdir, "sce_DLPFC_annotated")
)
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

# Gene filtering
genes_to_use <- rownames(spe)
genes_to_use <- setdiff(genes_to_use,
    grep("^MT-|^mt-", genes_to_use, value = TRUE))
genes_to_use <- setdiff(genes_to_use,
    grep("^RP[SL]|^Rp[sl]", genes_to_use, value = TRUE))
spe_sub <- spe[genes_to_use, ]

# Normalize
spe_sub <- scater::logNormCounts(spe_sub)

# HVGs
dec <- scran::modelGeneVar(spe_sub, assay.type = "logcounts")
hvgs <- scran::getTopHVGs(dec, n = 2000)

cat("  Matrix dimensions (all genes):", dim(spe_sub), "\n")
cat("  HVGs:", length(hvgs), "\n")

# =============================================================================
# 2. Check matrix class/sparsity
# =============================================================================

cat("\n=== Matrix diagnostics ===\n")
lc <- assay(spe_sub, "logcounts")
cat("  logcounts class:", class(lc), "\n")
cat("  logcounts dimensions:", dim(lc), "\n")

lc_hvg <- lc[hvgs, ]
cat("  logcounts[HVGs,] class:", class(lc_hvg), "\n")
cat("  logcounts[HVGs,] dimensions:", dim(lc_hvg), "\n")

if (is(lc_hvg, "sparseMatrix")) {
    nnz <- length(lc_hvg@x)
    total <- prod(dim(lc_hvg))
    cat("  Sparsity:", round(1 - nnz / total, 3), "\n")
    cat("  Non-zeros:", nnz, "/", total, "\n")
} else {
    cat("  Matrix is dense!\n")
}

# =============================================================================
# 3. Benchmark different PCA methods via scater::runPCA
# =============================================================================

cat("\n=== PCA benchmarks via scater::runPCA ===\n\n")

# Exact SVD (default)
cat("--- runPCA default (exact SVD) ---\n")
t1 <- system.time({
    res1 <- scater::runPCA(spe_sub, subset_row = hvgs,
                            ncomponents = 30, exprs_values = "logcounts")
})
cat("  Time:", t1["elapsed"], "s\n\n")

# IrlbaParam
cat("--- runPCA + IrlbaParam ---\n")
t2 <- system.time({
    res2 <- scater::runPCA(spe_sub, subset_row = hvgs,
                            ncomponents = 30, exprs_values = "logcounts",
                            BSPARAM = IrlbaParam())
})
cat("  Time:", t2["elapsed"], "s\n\n")

# RandomParam (randomized SVD — often fastest for this shape)
cat("--- runPCA + RandomParam ---\n")
t3 <- system.time({
    res3 <- scater::runPCA(spe_sub, subset_row = hvgs,
                            ncomponents = 30, exprs_values = "logcounts",
                            BSPARAM = RandomParam())
})
cat("  Time:", t3["elapsed"], "s\n\n")

# IrlbaParam with deferred centering (avoids densifying)
cat("--- runPCA + IrlbaParam(deferred=TRUE) ---\n")
t4 <- system.time({
    res4 <- scater::runPCA(spe_sub, subset_row = hvgs,
                            ncomponents = 30, exprs_values = "logcounts",
                            BSPARAM = IrlbaParam(deferred = TRUE))
})
cat("  Time:", t4["elapsed"], "s\n\n")

# RandomParam with deferred centering
cat("--- runPCA + RandomParam(deferred=TRUE) ---\n")
t5 <- system.time({
    res5 <- scater::runPCA(spe_sub, subset_row = hvgs,
                            ncomponents = 30, exprs_values = "logcounts",
                            BSPARAM = RandomParam(deferred = TRUE))
})
cat("  Time:", t5["elapsed"], "s\n\n")

# =============================================================================
# 4. Manual PCA to isolate centering vs SVD
# =============================================================================

cat("=== Manual PCA: isolating centering vs SVD ===\n\n")

mat_hvg <- lc_hvg  # sparse, 2000 x 77604

# Time: densify + center
cat("--- Realize to dense + center ---\n")
t_dense <- system.time({
    mat_dense <- as.matrix(Matrix::t(mat_hvg))  # 77604 x 2000
    col_means <- colMeans(mat_dense)
    mat_centered <- sweep(mat_dense, 2, col_means)
})
cat("  Time:", t_dense["elapsed"], "s\n")
cat("  Dense matrix size:", format(object.size(mat_centered), units = "GB"), "\n\n")

# Time: SVD on dense centered matrix
cat("--- Exact SVD (base::svd, 30 components) ---\n")
t_svd <- system.time({
    sv <- svd(mat_centered, nu = 30, nv = 30)
})
cat("  Time:", t_svd["elapsed"], "s\n\n")

# Time: irlba on dense centered matrix
cat("--- irlba::irlba on dense centered matrix ---\n")
t_irlba_dense <- system.time({
    ir <- irlba::irlba(mat_centered, nv = 30)
})
cat("  Time:", t_irlba_dense["elapsed"], "s\n\n")

# Time: irlba on sparse matrix with center argument
cat("--- irlba::irlba on sparse matrix (deferred center) ---\n")
mat_sparse_t <- Matrix::t(mat_hvg)  # 77604 x 2000
col_means_sparse <- Matrix::colMeans(mat_sparse_t)
t_irlba_sparse <- system.time({
    ir2 <- irlba::irlba(mat_sparse_t, nv = 30, center = col_means_sparse)
})
cat("  Time:", t_irlba_sparse["elapsed"], "s\n\n")

# =============================================================================
# 5. Summary
# =============================================================================

cat("====================================================\n")
cat("SUMMARY (seconds)\n")
cat("====================================================\n")
timings <- data.frame(
    Method = c(
        "runPCA default (exact)",
        "runPCA + IrlbaParam",
        "runPCA + RandomParam",
        "runPCA + IrlbaParam(deferred=TRUE)",
        "runPCA + RandomParam(deferred=TRUE)",
        "Manual: realize dense + center",
        "Manual: base::svd (dense)",
        "Manual: irlba (dense centered)",
        "Manual: irlba (sparse + deferred center)"
    ),
    Seconds = c(
        t1["elapsed"], t2["elapsed"], t3["elapsed"],
        t4["elapsed"], t5["elapsed"],
        t_dense["elapsed"], t_svd["elapsed"],
        t_irlba_dense["elapsed"], t_irlba_sparse["elapsed"]
    )
)
print(timings, row.names = FALSE)

# Check that deferred results are similar
cat("\n--- Correlation of PC1 across methods ---\n")
pc1_exact <- reducedDim(res1, "PCA")[, 1]
pc1_irlba <- reducedDim(res2, "PCA")[, 1]
pc1_random <- reducedDim(res3, "PCA")[, 1]
pc1_irlba_def <- reducedDim(res4, "PCA")[, 1]
pc1_random_def <- reducedDim(res5, "PCA")[, 1]

cat("  exact vs irlba:         ", round(abs(cor(pc1_exact, pc1_irlba)), 4), "\n")
cat("  exact vs random:        ", round(abs(cor(pc1_exact, pc1_random)), 4), "\n")
cat("  exact vs irlba(deferred):", round(abs(cor(pc1_exact, pc1_irlba_def)), 4), "\n")
cat("  exact vs random(deferred):", round(abs(cor(pc1_exact, pc1_random_def)), 4), "\n")

cat("\nPCA benchmark complete.\n")
