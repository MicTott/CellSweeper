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

    # Injected outliers have near-zero library size <U+2014> should be flagged
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

# ============================================================
# Performance parameter tests
# ============================================================

test_that("clusterLocalOutliers BNPARAM is accepted", {
    result <- clusterLocalOutliers(spe,
        metrics = c("sum"),
        morpho_metrics = NULL,
        k = 20,
        BNPARAM = BiocNeighbors::KmknnParam())
    expect_true("cellsweeper_outlier" %in% colnames(colData(result)))
})

# ============================================================
# Modified-z formula correctness (regression test for the
# constant=1.4826 / 0.6745 inconsistency in 0.99.0)
# ============================================================

test_that("robust_z='mad' matches (x - median) / mad(x) convention", {
    # Build a single-cluster SPE with a hand-computable z-score for one cell.
    n <- 60
    coords <- cbind(runif(n), runif(n))
    counts <- matrix(rpois(20 * n, 5), nrow = 20, ncol = n)
    rownames(counts) <- paste0("g", seq_len(20))
    colnames(counts) <- paste0("c", seq_len(n))
    s <- SpatialExperiment(
        assays        = list(counts = as(counts, "CsparseMatrix")),
        spatialCoords = coords,
        sample_id     = rep("s", n)
    )
    s$sum      <- c(rnorm(n - 1, mean = 100, sd = 5), 200) # last cell extreme
    s$detected <- s$sum
    s$cell_cluster <- rep("A", n)

    # k=n-1 so every cell sees every other cell as neighbour ->
    # neighbourhood = full cluster minus self.
    res <- clusterLocalOutliers(s, metrics = "sum", morpho_metrics = NULL,
        k = n - 1, log_transform = FALSE,
        robust_z = "mad", direction = "higher")

    other <- s$sum[-n]
    expected <- (s$sum[n] - stats::median(other)) / stats::mad(other)
    expect_equal(colData(res)$sum_zscore[n], expected, tolerance = 1e-9)
    expect_true(colData(res)$sum_outlier[n])
})

test_that("robust_z='iglewicz' matches 0.6745 * (x - median) / MAD_raw", {
    n <- 60
    coords <- cbind(runif(n), runif(n))
    counts <- matrix(rpois(20 * n, 5), nrow = 20, ncol = n)
    rownames(counts) <- paste0("g", seq_len(20))
    colnames(counts) <- paste0("c", seq_len(n))
    s <- SpatialExperiment(
        assays        = list(counts = as(counts, "CsparseMatrix")),
        spatialCoords = coords,
        sample_id     = rep("s", n)
    )
    s$sum <- c(rnorm(n - 1, mean = 100, sd = 5), 200)
    s$cell_cluster <- rep("A", n)

    res <- clusterLocalOutliers(s, metrics = "sum", morpho_metrics = NULL,
        k = n - 1, log_transform = FALSE,
        robust_z = "iglewicz", direction = "higher")
    other <- s$sum[-n]
    expected <- 0.6745 *
        (s$sum[n] - stats::median(other)) /
        stats::mad(other, constant = 1)
    expect_equal(colData(res)$sum_zscore[n], expected, tolerance = 1e-9)
})
