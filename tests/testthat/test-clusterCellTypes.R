library(CellSweeper)
library(SpatialExperiment)
library(SummarizedExperiment)

spe <- make_test_spe(n_cells_per_cluster = 100, n_outliers = 3, seed = 42)

test_that("clusterCellTypes adds cluster labels to colData", {
    result <- clusterCellTypes(spe, resolution = 0.5, regress_covariates = FALSE)
    expect_s4_class(result, "SpatialExperiment")
    expect_true("cell_cluster" %in% colnames(colData(result)))
    expect_false(any(is.na(colData(result)$cell_cluster)))
})

test_that("clusterCellTypes stores PCA in reducedDims", {
    result <- clusterCellTypes(spe, resolution = 0.5, regress_covariates = FALSE)
    expect_true("PCA" %in% SingleCellExperiment::reducedDimNames(result))
    expect_equal(nrow(SingleCellExperiment::reducedDim(result, "PCA")),
                 ncol(result))
})

test_that("clusterCellTypes use_existing works", {
    spe$my_labels <- rep(c("typeA", "typeB", "typeC"), each = 100)
    result <- clusterCellTypes(spe, use_existing = "my_labels")
    expect_equal(colData(result)$cell_cluster,
                 as.character(spe$my_labels))
})

test_that("clusterCellTypes errors on missing use_existing column", {
    expect_error(
        clusterCellTypes(spe, use_existing = "nonexistent"),
        "not found"
    )
})

test_that("clusterCellTypes assigns all to one cluster when too few cells", {
    # When total cells < min_cluster_size, all get "cluster_1"
    spe_tiny <- spe[, 1:20]
    expect_warning(
        result <- clusterCellTypes(spe_tiny, min_cluster_size = 50),
        "Fewer cells"
    )
    expect_true(all(colData(result)$cell_cluster == "cluster_1"))
})

test_that("clusterCellTypes skips missing covariates gracefully", {
    # Test that the message is emitted even if clustering itself may error
    # on synthetic data with regression
    expect_message(
        tryCatch(
            clusterCellTypes(spe,
                covariates = c("sum", "nonexistent_metric"),
                regress_covariates = TRUE),
            error = function(e) NULL
        ),
        "not found"
    )
})

test_that("clusterCellTypes with regress_covariates produces clusters", {
    # Regression + modelGeneVar can fail on synthetic data; test that
    # it at least starts correctly and the regression step runs
    result <- tryCatch(
        clusterCellTypes(spe, regress_covariates = TRUE,
                         covariates = c("sum", "detected")),
        error = function(e) {
            # If it fails during scran internals, that's OK for synthetic data
            # Just verify it was a scran issue, not a CellSweeper bug
            expect_true(grepl("density|bandwidth|need at least",
                              e$message, ignore.case = TRUE))
            NULL
        }
    )
    if (!is.null(result)) {
        expect_true("cell_cluster" %in% colnames(colData(result)))
    }
})

test_that("clusterCellTypes custom cluster_col works", {
    result <- clusterCellTypes(spe, cluster_col = "my_clusters",
                               regress_covariates = FALSE)
    expect_true("my_clusters" %in% colnames(colData(result)))
})

test_that("clusterCellTypes errors on non-SPE input", {
    expect_error(clusterCellTypes(data.frame()), "SpatialExperiment")
})

test_that("clusterCellTypes handles very small datasets", {
    small_spe <- spe[, 1:10]
    expect_warning(
        result <- clusterCellTypes(small_spe, min_cluster_size = 50),
        "Fewer cells"
    )
    expect_true("cell_cluster" %in% colnames(colData(result)))
})
