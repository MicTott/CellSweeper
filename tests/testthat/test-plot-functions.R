library(CellSweeper)
library(SpatialExperiment)
library(SummarizedExperiment)

spe <- make_test_spe(n_cells_per_cluster = 50, n_outliers = 3, seed = 42)
spe$cell_cluster <- spe$true_cluster

# Add mock QC columns that plotting functions expect
colData(spe)$sum_zscore <- rnorm(ncol(spe))
colData(spe)$sum_outlier <- abs(colData(spe)$sum_zscore) > 2
colData(spe)$cellsweeper_outlier <- colData(spe)$sum_outlier
colData(spe)$cluster_qc_mahal <- runif(ncol(spe), 0, 10)
colData(spe)$cluster_spatial_homogeneity <- runif(ncol(spe), 0.5, 1)
colData(spe)$cluster_artifact_flag <- FALSE

# ============================================================
# Point-based tests (original)
# ============================================================

test_that("plotSpatialOutliers returns ggplot", {
    skip_if_not_installed("ggplot2")
    p <- plotSpatialOutliers(spe, metric = "sum")
    expect_s3_class(p, "ggplot")
})

test_that("plotClusterQC returns ggplot", {
    skip_if_not_installed("ggplot2")
    p <- plotClusterQC(spe, metric = "sum")
    expect_s3_class(p, "ggplot")
})

test_that("plotSpatialDispersion returns ggplot", {
    skip_if_not_installed("ggplot2")
    p <- plotSpatialDispersion(spe)
    expect_s3_class(p, "ggplot")
})

test_that("plotClusterSummary returns list of ggplots", {
    skip_if_not_installed("ggplot2")
    result <- plotClusterSummary(spe, metrics = c("sum"))
    expect_type(result, "list")
    if ("mahalanobis" %in% names(result)) {
        expect_s3_class(result$mahalanobis, "ggplot")
    }
    if ("homogeneity" %in% names(result)) {
        expect_s3_class(result$homogeneity, "ggplot")
    }
})

test_that("plotSpatialOutliers errors on missing metric", {
    skip_if_not_installed("ggplot2")
    expect_error(plotSpatialOutliers(spe, metric = "nonexistent"),
                 "not found")
})

test_that("plotClusterQC errors on missing cluster column", {
    skip_if_not_installed("ggplot2")
    expect_error(plotClusterQC(spe, cluster_col = "nonexistent"),
                 "not found")
})

test_that("plotSpatialOutliers errors on non-SPE", {
    skip_if_not_installed("ggplot2")
    expect_error(plotSpatialOutliers(data.frame()), "SpatialExperiment")
})

# ============================================================
# Polygon auto-detection and fallback tests
# ============================================================

test_that("plotSpatialOutliers auto-falls back to points when no polygons", {
    skip_if_not_installed("ggplot2")
    p <- plotSpatialOutliers(spe, metric = "sum", use_polygons = "auto")
    expect_s3_class(p, "ggplot")
    layer_geoms <- vapply(p$layers, function(l) class(l$geom)[1],
                          character(1))
    expect_true("GeomPoint" %in% layer_geoms)
})

test_that("plotSpatialOutliers errors when use_polygons=TRUE but none exist", {
    skip_if_not_installed("ggplot2")
    expect_error(
        plotSpatialOutliers(spe, metric = "sum", use_polygons = TRUE),
        "not found")
})

test_that("plotSpatialOutliers use_polygons=FALSE forces points", {
    skip_if_not_installed("ggplot2")
    p <- plotSpatialOutliers(spe, metric = "sum", use_polygons = FALSE)
    expect_s3_class(p, "ggplot")
})

test_that("plotSpatialDispersion auto-falls back to points when no polygons", {
    skip_if_not_installed("ggplot2")
    p <- plotSpatialDispersion(spe, use_polygons = "auto")
    expect_s3_class(p, "ggplot")
})

test_that("plotSpatialDispersion errors when use_polygons=TRUE but none", {
    skip_if_not_installed("ggplot2")
    expect_error(
        plotSpatialDispersion(spe, use_polygons = TRUE),
        "not found")
})

# ============================================================
# Polygon rendering tests (require sf)
# ============================================================

test_that("plotSpatialOutliers renders polygons when available", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("sf")

    spe_poly <- make_test_spe(n_cells_per_cluster = 50, n_outliers = 3,
                              seed = 42, include_polygons = TRUE)
    colData(spe_poly)$cellsweeper_outlier <-
        abs(rnorm(ncol(spe_poly))) > 2

    p <- plotSpatialOutliers(spe_poly, metric = "sum",
                             use_polygons = TRUE)
    expect_s3_class(p, "ggplot")
    layer_geoms <- vapply(p$layers, function(l) class(l$geom)[1],
                          character(1))
    expect_true("GeomSf" %in% layer_geoms)
})

test_that("plotSpatialDispersion renders polygons when available", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("sf")

    spe_poly <- make_test_spe(n_cells_per_cluster = 50, n_outliers = 3,
                              seed = 42, include_polygons = TRUE)
    spe_poly$cell_cluster <- spe_poly$true_cluster

    p <- plotSpatialDispersion(spe_poly, use_polygons = TRUE)
    expect_s3_class(p, "ggplot")
    layer_geoms <- vapply(p$layers, function(l) class(l$geom)[1],
                          character(1))
    expect_true("GeomSf" %in% layer_geoms)
})

test_that(".hasPolygons returns FALSE for missing column", {
    expect_false(CellSweeper:::.hasPolygons(spe, "nonexistent"))
})

test_that(".extractPolygonData returns NULL for missing column", {
    result <- CellSweeper:::.extractPolygonData(spe, "nonexistent", "sum")
    expect_null(result)
})
