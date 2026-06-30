library(CellSweeper)
library(SpatialExperiment)
library(SummarizedExperiment)

spe <- make_test_spe(n_cells_per_cluster = 100, n_outliers = 5, seed = 42)

test_that("runCellSweeper full pipeline produces expected columns", {
    # Use pre-existing labels to skip clustering (faster test)
    result <- runCellSweeper(spe,
        use_existing = "true_cluster",
        metrics = c("sum", "detected"),
        morpho_metrics = NULL,
        k = 20,
        n_permutations = 20,
        verbose = FALSE)

    expect_s4_class(result, "SpatialExperiment")

    # Check all expected output columns
    expected <- c("cell_cluster",
        "cluster_qc_mahal", "cluster_spatial_homogeneity",
        "cluster_spatial_pvalue", "cluster_artifact_flag",
        "sum_zscore", "sum_outlier",
        "detected_zscore", "detected_outlier",
        "cellsweeper_outlier", "cellsweeper_qc")

    for (col in expected) {
        expect_true(col %in% colnames(colData(result)),
            info = paste("Missing column:", col))
    }
})

test_that("cellsweeper_qc is union of artifact and outlier flags", {
    result <- runCellSweeper(spe,
        use_existing = "true_cluster",
        metrics = c("sum"),
        morpho_metrics = NULL,
        k = 20,
        n_permutations = 20,
        verbose = FALSE)

    cd <- colData(result)
    # cellsweeper_qc should be TRUE when either flag is TRUE
    expected_qc <- cd$cluster_artifact_flag | cd$cellsweeper_outlier
    expected_qc[is.na(expected_qc)] <- FALSE
    expect_equal(cd$cellsweeper_qc, expected_qc)
})

test_that("runCellSweeper removes global filter cells", {
    # With strict thresholds, some cells should be removed
    result <- runCellSweeper(spe,
        min_counts = 10,
        use_existing = "true_cluster",
        metrics = c("sum"),
        morpho_metrics = NULL,
        k = 20,
        n_permutations = 20,
        verbose = FALSE)

    expect_lt(ncol(result), ncol(spe))
})

test_that("runCellSweeper errors on non-SPE input", {
    expect_error(runCellSweeper(data.frame()), "SpatialExperiment")
})
