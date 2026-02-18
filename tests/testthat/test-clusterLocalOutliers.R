library(CellSweeper)
library(SpatialExperiment)
library(SummarizedExperiment)

spe <- make_test_spe(n_cells_per_cluster = 100, n_outliers = 5, seed = 42)
spe$cell_cluster <- spe$true_cluster

test_that("clusterLocalOutliers adds expected colData columns", {
    result <- clusterLocalOutliers(spe,
                                    metrics = c("sum", "detected"),
                                    morpho_metrics = NULL,
                                    k = 20)

    expect_true("sum_zscore" %in% colnames(colData(result)))
    expect_true("sum_outlier" %in% colnames(colData(result)))
    expect_true("detected_zscore" %in% colnames(colData(result)))
    expect_true("detected_outlier" %in% colnames(colData(result)))
    expect_true("cellsweeper_outlier" %in% colnames(colData(result)))
})

test_that("clusterLocalOutliers detects injected outliers", {
    result <- clusterLocalOutliers(spe,
                                    metrics = c("sum"),
                                    morpho_metrics = NULL,
                                    k = 20, z_threshold = 3)

    # Injected outliers have near-zero library size — should be flagged
    flagged <- colData(result)$cellsweeper_outlier
    expect_true(any(flagged, na.rm = TRUE))
    # Should not flag most cells
    expect_lt(sum(flagged, na.rm = TRUE), ncol(spe) / 2)
})

test_that("clusterLocalOutliers respects direction parameter", {
    # Lower direction: flags cells with LOW values
    result_lower <- clusterLocalOutliers(spe,
                                          metrics = c("sum"),
                                          morpho_metrics = NULL,
                                          direction = "lower",
                                          k = 20)

    # Higher direction: flags cells with HIGH values
    result_higher <- clusterLocalOutliers(spe,
                                           metrics = c("sum"),
                                           morpho_metrics = NULL,
                                           direction = "higher",
                                           k = 20)

    # With our injected low-count outliers, "lower" should flag more
    n_lower <- sum(colData(result_lower)$cellsweeper_outlier, na.rm = TRUE)
    n_higher <- sum(colData(result_higher)$cellsweeper_outlier, na.rm = TRUE)
    expect_gte(n_lower, n_higher)
})

test_that("union combines more aggressively than intersection", {
    result_union <- clusterLocalOutliers(spe,
                                          metrics = c("sum", "detected"),
                                          morpho_metrics = NULL,
                                          combine = "union", k = 20)
    result_inter <- clusterLocalOutliers(spe,
                                          metrics = c("sum", "detected"),
                                          morpho_metrics = NULL,
                                          combine = "intersection", k = 20)

    n_union <- sum(colData(result_union)$cellsweeper_outlier, na.rm = TRUE)
    n_inter <- sum(colData(result_inter)$cellsweeper_outlier, na.rm = TRUE)
    expect_gte(n_union, n_inter)
})

test_that("clusterLocalOutliers handles morpho_metrics gracefully", {
    # With morpho metrics present
    result <- clusterLocalOutliers(spe,
                                    metrics = c("sum"),
                                    morpho_metrics = c("Area_um"),
                                    k = 20)
    expect_true("Area_um_zscore" %in% colnames(colData(result)))

    # With morpho metrics absent
    spe_no_morpho <- spe
    colData(spe_no_morpho)$Area_um <- NULL
    expect_message(
        result2 <- clusterLocalOutliers(spe_no_morpho,
                                         metrics = c("sum"),
                                         morpho_metrics = c("Area_um"),
                                         k = 20),
        "not found"
    )
    expect_true("cellsweeper_outlier" %in% colnames(colData(result2)))
})

test_that("clusterLocalOutliers skips artifact clusters", {
    spe_art <- spe
    # Mark cluster_1 as artifact
    spe_art$cluster_artifact_flag <-
        spe_art$cell_cluster == "cluster_1"

    result <- clusterLocalOutliers(spe_art,
                                    metrics = c("sum"),
                                    morpho_metrics = NULL,
                                    exclude_artifact_clusters = TRUE,
                                    k = 20)

    # Cells in artifact cluster should have NA z-scores
    art_zscores <- colData(result)$sum_zscore[
        spe_art$cell_cluster == "cluster_1"]
    expect_true(all(is.na(art_zscores)))
})

test_that("clusterLocalOutliers errors on missing cluster column", {
    expect_error(
        clusterLocalOutliers(spe, cluster_col = "nonexistent"),
        "not found"
    )
})

test_that("clusterLocalOutliers handles no sample_id", {
    spe_no_sample <- spe
    colData(spe_no_sample)$sample_id <- NULL
    result <- clusterLocalOutliers(spe_no_sample,
                                    metrics = c("sum"),
                                    morpho_metrics = NULL,
                                    k = 20)
    expect_true("cellsweeper_outlier" %in% colnames(colData(result)))
})

test_that("z-scores are numeric and mostly finite", {
    result <- clusterLocalOutliers(spe,
                                    metrics = c("sum"),
                                    morpho_metrics = NULL,
                                    k = 20)
    z <- colData(result)$sum_zscore
    expect_true(is.numeric(z))
    # Most should be finite (NAs only for cells in skipped clusters)
    n_finite <- sum(is.finite(z), na.rm = TRUE)
    expect_gt(n_finite, ncol(spe) / 2)
})
