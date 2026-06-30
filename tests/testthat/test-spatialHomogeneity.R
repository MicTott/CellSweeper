library(CellSweeper)
library(SpatialExperiment)
library(SummarizedExperiment)

spe <- make_test_spe(n_cells_per_cluster = 100, n_outliers = 0, seed = 42)
# Use known spatial clusters as cluster labels
spe$cell_cluster <- spe$true_cluster

test_that("computeNeighborhoodHomogeneity adds column to colData", {
    result <- computeNeighborhoodHomogeneity(spe)
    expect_true("neighborhood_homogeneity" %in% colnames(colData(result)))
    expect_equal(length(colData(result)$neighborhood_homogeneity), ncol(result))
})

test_that("spatially coherent clusters have high homogeneity", {
    # Our synthetic data has 3 well-separated spatial blobs
    result <- computeNeighborhoodHomogeneity(spe, k = 20)
    homo <- colData(result)$neighborhood_homogeneity

    # Mean homogeneity should be high since clusters are spatially separated
    expect_gt(mean(homo, na.rm = TRUE), 0.5)
})

test_that("random labels have lower homogeneity than real clusters", {
    # Real cluster labels
    result_real <- computeNeighborhoodHomogeneity(spe, k = 20)
    homo_real <- mean(colData(result_real)$neighborhood_homogeneity,
        na.rm = TRUE)

    # Random labels
    spe_random <- spe
    spe_random$cell_cluster <- sample(spe$cell_cluster)
    result_random <- computeNeighborhoodHomogeneity(spe_random, k = 20)
    homo_random <- mean(colData(result_random)$neighborhood_homogeneity,
        na.rm = TRUE)

    expect_gt(homo_real, homo_random)
})

test_that("computeNeighborhoodHomogeneity errors on missing cluster_col", {
    expect_error(
        computeNeighborhoodHomogeneity(spe, cluster_col = "nonexistent"),
        "not found"
    )
})

test_that("permuteHomogeneity returns correct structure", {
    perm_results <- permuteHomogeneity(spe, k = 10, n_permutations = 20)
    expect_s3_class(perm_results, "data.frame")
    expect_true(all(c("sample", "cluster", "observed_homogeneity",
        "null_mean", "null_sd", "z_score", "p_value") %in%
        colnames(perm_results)))
    # One row per (sample, cluster) pair <U+2014> single sample, so same as n clusters
    expect_equal(nrow(perm_results),
        length(unique(spe$cell_cluster)))
})

test_that("spatially coherent clusters have low p-values", {
    perm_results <- permuteHomogeneity(spe, k = 10, n_permutations = 100)
    # Well-separated spatial blobs should have very low p-values
    expect_true(all(perm_results$p_value < 0.1))
})

test_that("permuteHomogeneity handles no sample_id column", {
    spe_no_sample <- spe
    colData(spe_no_sample)$sample_id <- NULL
    # Should work by treating all cells as one sample
    perm_results <- permuteHomogeneity(spe_no_sample, k = 10,
        n_permutations = 10)
    expect_s3_class(perm_results, "data.frame")
})

# ============================================================
# Performance parameter tests
# ============================================================

test_that("computeNeighborhoodHomogeneity accepts BNPARAM", {
    result <- computeNeighborhoodHomogeneity(spe, k = 10,
        BNPARAM = BiocNeighbors::KmknnParam())
    expect_true("neighborhood_homogeneity" %in% colnames(colData(result)))
})

test_that("permuteHomogeneity accepts BNPARAM", {
    perm_results <- permuteHomogeneity(spe, k = 10, n_permutations = 5,
        BNPARAM = BiocNeighbors::KmknnParam())
    expect_s3_class(perm_results, "data.frame")
    expect_true("p_value" %in% colnames(perm_results))
})

test_that("permuteHomogeneity accepts BPPARAM with SerialParam", {
    skip_if_not_installed("BiocParallel")
    perm_results <- permuteHomogeneity(spe, k = 10, n_permutations = 5,
        BPPARAM = BiocParallel::SerialParam())
    expect_s3_class(perm_results, "data.frame")
    expect_equal(nrow(perm_results), length(unique(spe$cell_cluster)))
})

test_that("permuteHomogeneity is reproducible across repeated runs", {
    skip_if_not_installed("BiocParallel")
    r1 <- permuteHomogeneity(spe, k = 10, n_permutations = 25,
        BPPARAM = BiocParallel::SerialParam(), seed = 7)
    r2 <- permuteHomogeneity(spe, k = 10, n_permutations = 25,
        BPPARAM = BiocParallel::SerialParam(), seed = 7)
    expect_equal(r1, r2)
})

test_that("permuteHomogeneity is reproducible regardless of BPPARAM", {
    skip_if_not_installed("BiocParallel")
    r_serial <- permuteHomogeneity(spe, k = 10, n_permutations = 25,
        BPPARAM = BiocParallel::SerialParam(), seed = 13)
    r_default <- permuteHomogeneity(spe, k = 10, n_permutations = 25,
        BPPARAM = NULL, seed = 13)
    expect_equal(r_serial$observed_homogeneity,
        r_default$observed_homogeneity)
    expect_equal(r_serial$p_value, r_default$p_value)
})

# ============================================================
# Multi-sample tests
# ============================================================

test_that("permuteHomogeneity returns per-(sample, cluster) rows", {
    spe_ms <- make_test_spe_multisample(n_cells_per_cluster = 50,
        n_outliers = 0)
    spe_ms$cell_cluster <- spe_ms$true_cluster
    perm_results <- permuteHomogeneity(spe_ms, k = 10, n_permutations = 10)

    expect_true("sample" %in% colnames(perm_results))
    # 2 samples x 3 clusters each = 6 rows
    n_sc <- length(unique(paste0(spe_ms$sample_id, ":", spe_ms$cell_cluster)))
    expect_equal(nrow(perm_results), n_sc)
})

test_that("computeNeighborhoodHomogeneity works with multi-sample data", {
    spe_ms <- make_test_spe_multisample(n_cells_per_cluster = 50,
        n_outliers = 0)
    spe_ms$cell_cluster <- spe_ms$true_cluster
    result <- computeNeighborhoodHomogeneity(spe_ms, k = 10)
    expect_true("neighborhood_homogeneity" %in% colnames(colData(result)))
    expect_equal(length(colData(result)$neighborhood_homogeneity), ncol(result))
    expect_false(any(is.na(colData(result)$neighborhood_homogeneity)))
})
