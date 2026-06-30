library(CellSweeper)
library(SpatialExperiment)
library(SummarizedExperiment)

spe <- make_test_spe(n_cells_per_cluster = 50, n_outliers = 3, seed = 42)
n_start <- ncol(spe)

test_that("globalFilter removes low-count cells", {
    # Outlier cells have ~2 counts, min_counts = 5 should remove them
    result <- globalFilter(spe, min_counts = 5, min_genes = 1)
    expect_s4_class(result, "SpatialExperiment")
    expect_lt(ncol(result), n_start)
    # All remaining cells should have >= 5 counts
    expect_true(all(result$sum >= 5))
})

test_that("globalFilter removes low-gene cells", {
    result <- globalFilter(spe, min_counts = 1, min_genes = 5)
    expect_s4_class(result, "SpatialExperiment")
    expect_true(all(result$detected >= 5))
})

test_that("globalFilter applies area filters when column exists", {
    result <- globalFilter(spe, min_counts = 1, min_genes = 1,
        min_area = 20, max_area = 200)
    expect_true(all(colData(result)$Area_um >= 20))
    expect_true(all(colData(result)$Area_um <= 200))
})

test_that("globalFilter skips area filter when column missing", {
    # Remove area column
    spe_no_area <- spe
    colData(spe_no_area)$Area_um <- NULL
    # Should still work, just skip area filtering
    expect_message(
        result <- globalFilter(spe_no_area, min_counts = 1, min_genes = 1,
            min_area = 20),
        "not found"
    )
    expect_s4_class(result, "SpatialExperiment")
})

test_that("globalFilter adds CountArea and log2CountArea when area exists", {
    # Remove pre-computed density columns so globalFilter computes them
    spe_tmp <- spe
    colData(spe_tmp)$CountArea <- NULL
    colData(spe_tmp)$log2CountArea <- NULL
    result <- globalFilter(spe_tmp, min_counts = 1, min_genes = 1)
    expect_true("CountArea" %in% colnames(colData(result)))
    expect_true("log2CountArea" %in% colnames(colData(result)))
    # CountArea should be counts / area
    expected <- result$sum / colData(result)$Area_um
    expect_equal(colData(result)$CountArea, expected)
    # log2CountArea should be log2(density + 1)
    expect_equal(colData(result)$log2CountArea, log2(expected + 1))
})

test_that("globalFilter skips density when log2CountArea already exists", {
    # helper-synthetic.R already adds CountArea + log2CountArea
    result <- globalFilter(spe, min_counts = 1, min_genes = 1)
    expect_true("log2CountArea" %in% colnames(colData(result)))
})

test_that("globalFilter does not add density when area column missing", {
    spe_no_area <- spe
    colData(spe_no_area)$Area_um <- NULL
    colData(spe_no_area)$CountArea <- NULL
    colData(spe_no_area)$log2CountArea <- NULL
    result <- globalFilter(spe_no_area, min_counts = 1, min_genes = 1)
    expect_false("CountArea" %in% colnames(colData(result)))
    expect_false("log2CountArea" %in% colnames(colData(result)))
})

test_that("globalFilter passes all cells when thresholds are permissive", {
    result <- globalFilter(spe, min_counts = 0, min_genes = 0)
    expect_equal(ncol(result), n_start)
})

test_that("globalFilter errors on non-SPE input", {
    expect_error(globalFilter(data.frame()), "SpatialExperiment")
})

test_that("globalFilter errors when QC columns missing", {
    spe_bare <- SpatialExperiment(
        assays = list(counts = matrix(1, nrow = 5, ncol = 5)),
        spatialCoords = matrix(rnorm(10), ncol = 2)
    )
    expect_error(globalFilter(spe_bare), "not found")
})
