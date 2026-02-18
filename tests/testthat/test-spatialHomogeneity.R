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
    expect_true(all(c("cluster", "observed_homogeneity", "null_mean",
                       "null_sd", "z_score", "p_value") %in%
                        colnames(perm_results)))
    # One row per unique cluster
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
