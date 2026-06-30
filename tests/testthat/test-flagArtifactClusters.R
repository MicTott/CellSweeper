library(CellSweeper)
library(SpatialExperiment)
library(SummarizedExperiment)

spe <- make_test_spe(n_cells_per_cluster = 100, n_outliers = 3, seed = 42)
spe$cell_cluster <- spe$true_cluster

test_that("flagArtifactClusters adds expected colData columns", {
    result <- flagArtifactClusters(spe,
        metrics = c("sum", "detected"),
        n_permutations = 20)

    expected_cols <- c("cluster_qc_mahal", "cluster_spatial_homogeneity",
        "cluster_spatial_pvalue", "cluster_artifact_flag")
    for (col in expected_cols) {
        expect_true(col %in% colnames(colData(result)),
            info = paste("Missing column:", col))
    }
})

test_that("flagArtifactClusters artifact_flag is logical", {
    result <- flagArtifactClusters(spe,
        metrics = c("sum", "detected"),
        n_permutations = 20)
    expect_type(colData(result)$cluster_artifact_flag, "logical")
})

test_that("spatially coherent clusters are not flagged", {
    # Our synthetic clusters are well-separated spatial blobs
    # They should NOT be flagged as artifacts
    result <- flagArtifactClusters(spe,
        metrics = c("sum", "detected"),
        n_permutations = 50,
        require_both = TRUE)
    # With require_both = TRUE, coherent clusters should not be flagged
    # even if they have different QC profiles
    flagged <- unique(colData(result)$cluster_artifact_flag[
        colData(result)$cluster_spatial_pvalue < 0.05])
    expect_false(any(flagged))
})

test_that("flagArtifactClusters skips missing metrics gracefully", {
    expect_message(
        result <- flagArtifactClusters(spe,
            metrics = c("sum", "nonexistent_metric"),
            n_permutations = 10),
        "not found"
    )
    expect_true("cluster_artifact_flag" %in% colnames(colData(result)))
})

test_that("flagArtifactClusters errors on missing cluster column", {
    expect_error(
        flagArtifactClusters(spe, cluster_col = "nonexistent"),
        "not found"
    )
})

test_that("flagArtifactClusters warns with < 3 clusters", {
    spe_2cl <- spe
    spe_2cl$cell_cluster <- ifelse(seq_len(ncol(spe)) <= 150, "A", "B")
    expect_warning(
        result <- flagArtifactClusters(spe_2cl,
            metrics = c("sum", "detected"),
            n_permutations = 10),
        "Fewer than 3"
    )
    # Should still complete (spatial analysis still runs)
    expect_true("cluster_artifact_flag" %in% colnames(colData(result)))
})

test_that("flagArtifactClusters errors when no valid metrics", {
    expect_error(
        flagArtifactClusters(spe, metrics = c("fake1", "fake2"),
            n_permutations = 10),
        "No valid metrics"
    )
})

test_that("require_both = FALSE is more aggressive", {
    result_both <- flagArtifactClusters(spe,
        metrics = c("sum", "detected"),
        n_permutations = 20,
        require_both = TRUE)
    result_either <- flagArtifactClusters(spe,
        metrics = c("sum", "detected"),
        n_permutations = 20,
        require_both = FALSE)
    # Either condition should flag >= as many cells as both conditions
    n_flagged_both <- sum(colData(result_both)$cluster_artifact_flag)
    n_flagged_either <- sum(colData(result_either)$cluster_artifact_flag)
    expect_gte(n_flagged_either, n_flagged_both)
})

# ============================================================
# Multi-sample tests
# ============================================================

test_that("flagArtifactClusters works with multi-sample data", {
    spe_ms <- make_test_spe_multisample(n_cells_per_cluster = 100,
        n_outliers = 3)
    spe_ms$cell_cluster <- spe_ms$true_cluster
    result <- flagArtifactClusters(spe_ms,
        metrics = c("sum", "detected"),
        n_permutations = 10)

    expected_cols <- c("cluster_qc_mahal", "cluster_spatial_homogeneity",
        "cluster_spatial_pvalue", "cluster_artifact_flag")
    for (col in expected_cols) {
        expect_true(col %in% colnames(colData(result)),
            info = paste("Missing column:", col))
    }
    expect_type(colData(result)$cluster_artifact_flag, "logical")
    # All cells should have values
    expect_equal(length(colData(result)$cluster_artifact_flag), ncol(spe_ms))
})

test_that("flagArtifactClusters per-sample Mahalanobis computes independently", {
    spe_ms <- make_test_spe_multisample(n_cells_per_cluster = 100,
        n_outliers = 3)
    spe_ms$cell_cluster <- spe_ms$true_cluster
    result <- flagArtifactClusters(spe_ms,
        metrics = c("sum", "detected"),
        n_permutations = 10)

    # Each sample's cells should have Mahalanobis values
    cd <- colData(result)
    s1_mahal <- cd$cluster_qc_mahal[cd$sample_id == "sample1"]
    s2_mahal <- cd$cluster_qc_mahal[cd$sample_id == "sample2"]
    # Both should have non-NA values (3 clusters each -> Mahalanobis runs)
    expect_false(all(is.na(s1_mahal)))
    expect_false(all(is.na(s2_mahal)))
})
