# Suppress R CMD check NOTEs for ggplot2 aes() variables
utils::globalVariables(c(".data", "display_cluster", "cluster_qc_mahal",
                          "cluster_artifact_flag",
                          "cluster_spatial_homogeneity"))

# ============================================================
# Internal polygon helpers
# ============================================================

#' Check whether polygon plotting is possible
#'
#' @param spe A SpatialExperiment object
#' @param polygon_col Name of the colData column with sf geometry
#' @return Logical
#' @noRd
.hasPolygons <- function(spe, polygon_col) {
    requireNamespace("sf", quietly = TRUE) &&
        polygon_col %in% colnames(colData(spe)) &&
        inherits(colData(spe)[[polygon_col]], c("sfc", "sf"))
}

#' Extract polygon geometry from colData for ggplot2
#'
#' Builds an sf data.frame from the polygon geometry stored in colData,
#' suitable for use with \code{ggplot2::geom_sf()}.
#'
#' @param spe A SpatialExperiment object
#' @param polygon_col Name of the colData column containing sf geometry
#' @param extra_cols Character vector of colData column names to include
#' @return An sf data.frame with geometry and requested columns, or NULL
#'   if extraction fails
#' @noRd
.extractPolygonData <- function(spe, polygon_col, extra_cols = character(0)) {
    if (!requireNamespace("sf", quietly = TRUE)) return(NULL)

    cd <- colData(spe)
    if (!polygon_col %in% colnames(cd)) return(NULL)

    poly_geom <- cd[[polygon_col]]
    if (!inherits(poly_geom, c("sfc", "sf"))) return(NULL)

    # Build a plain data.frame with the requested columns
    available <- intersect(extra_cols, colnames(cd))
    extra_data <- as.data.frame(cd[, available, drop = FALSE])

    # Construct sf data.frame
    sf::st_sf(extra_data, geometry = poly_geom)
}


#' Plot Spatial Outliers
#'
#' Creates a spatial plot colored by a QC metric, with outlier cells
#' highlighted. Renders cell boundaries as filled polygons when available
#' (e.g., from SpaceTrooper preprocessing), otherwise falls back to points.
#'
#' @param spe A \linkS4class{SpatialExperiment} object with QC results from
#'   \code{\link{clusterLocalOutliers}}.
#' @param metric Name of the QC metric to color by (default "sum").
#' @param outlier_col Name of the \code{colData} column with outlier flags
#'   (default "cellsweeper_outlier").
#' @param sample Character. Which sample to plot (default NULL plots first
#'   sample).
#' @param samples Column name for sample IDs (default "sample_id").
#' @param point_size Point size for point-based rendering (default 0.5).
#' @param colors Color gradient for the metric (default
#'   \code{c("grey90", "navy")}).
#' @param outlier_color Color for outlier highlighting (default "red").
#' @param outlier_stroke Stroke width for outlier highlighting (default 0.5).
#' @param use_polygons Controls polygon rendering. \code{"auto"} (default)
#'   uses polygons if available in \code{colData} and the \pkg{sf} package is
#'   installed, otherwise falls back to points. \code{TRUE} requires polygons
#'   (errors if unavailable). \code{FALSE} always uses points.
#' @param polygon_col Name of the \code{colData} column containing \pkg{sf}
#'   polygon geometry (default \code{"polygons"}, matching SpaceTrooper
#'   output).
#' @param polygon_border_color Border color for polygons (default \code{NA},
#'   no borders).
#' @param polygon_border_width Border line width for polygons (default 0.1).
#'
#' @return A ggplot object.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' library(SpatialExperiment)
#' library(scuttle)
#'
#' counts <- matrix(rpois(2000, 5), nrow = 200, ncol = 10)
#' rownames(counts) <- paste0("Gene", 1:200)
#' coords <- matrix(rnorm(20), ncol = 2)
#' colnames(coords) <- c("x", "y")
#' spe <- SpatialExperiment(assays = list(counts = counts),
#'     spatialCoords = coords)
#' spe <- addPerCellQCMetrics(spe)
#'
#' plotSpatialOutliers(spe, metric = "sum")
plotSpatialOutliers <- function(spe,
                                metric = "sum",
                                outlier_col = "cellsweeper_outlier",
                                sample = NULL,
                                samples = "sample_id",
                                point_size = 0.5,
                                colors = c("grey90", "navy"),
                                outlier_color = "red",
                                outlier_stroke = 0.5,
                                use_polygons = "auto",
                                polygon_col = "polygons",
                                polygon_border_color = NA,
                                polygon_border_width = 0.1) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting. ",
             "Install with: install.packages('ggplot2')")
    }

    if (!is(spe, "SpatialExperiment")) {
        stop("'spe' must be a SpatialExperiment object.")
    }

    # --- Resolve polygon vs point mode ---
    render_polygons <- FALSE
    if (identical(use_polygons, TRUE)) {
        if (!.hasPolygons(spe, polygon_col)) {
            stop("Polygon column '", polygon_col,
                 "' not found or sf package not available.")
        }
        render_polygons <- TRUE
    } else if (identical(use_polygons, "auto")) {
        render_polygons <- .hasPolygons(spe, polygon_col)
    }

    cd <- as.data.frame(colData(spe))

    if (!metric %in% colnames(cd)) {
        stop("Metric '", metric, "' not found in colData.")
    }

    has_outliers <- outlier_col %in% colnames(cd) &&
        any(cd[[outlier_col]], na.rm = TRUE)

    if (render_polygons) {
        # --- Polygon-based rendering ---
        needed <- c(metric, samples)
        if (outlier_col %in% colnames(cd)) needed <- c(needed, outlier_col)
        sf_df <- .extractPolygonData(spe, polygon_col, needed)
        if (is.null(sf_df)) {
            render_polygons <- FALSE
        }
    }

    if (render_polygons) {
        # Sample subsetting
        if (samples %in% colnames(sf_df)) {
            if (is.null(sample)) sample <- sf_df[[samples]][1]
            sf_df <- sf_df[sf_df[[samples]] == sample, ]
        }

        p <- ggplot2::ggplot(sf_df) +
            ggplot2::geom_sf(
                ggplot2::aes(fill = .data[[metric]]),
                color = polygon_border_color,
                linewidth = polygon_border_width) +
            ggplot2::scale_fill_gradient(low = colors[1], high = colors[2]) +
            ggplot2::coord_sf() +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = paste0("Spatial: ", metric),
                          fill = metric)

        if (has_outliers) {
            sf_outlier <- sf_df[sf_df[[outlier_col]] == TRUE &
                                !is.na(sf_df[[outlier_col]]), ]
            if (nrow(sf_outlier) > 0) {
                p <- p + ggplot2::geom_sf(
                    data = sf_outlier,
                    fill = NA,
                    color = outlier_color,
                    linewidth = outlier_stroke)
            }
        }
    } else {
        # --- Point-based rendering (original) ---
        coords <- as.data.frame(spatialCoords(spe))
        df <- cbind(coords, cd)

        if (samples %in% colnames(df)) {
            if (is.null(sample)) sample <- df[[samples]][1]
            df <- df[df[[samples]] == sample, ]
        }

        coord_cols <- colnames(coords)

        p <- ggplot2::ggplot(df, ggplot2::aes(
                x = .data[[coord_cols[1]]],
                y = .data[[coord_cols[2]]])) +
            ggplot2::geom_point(
                ggplot2::aes(color = .data[[metric]]),
                size = point_size) +
            ggplot2::scale_color_gradient(low = colors[1],
                                          high = colors[2]) +
            ggplot2::coord_fixed() +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = paste0("Spatial: ", metric),
                          color = metric)

        if (has_outliers) {
            df_outlier <- df[df[[outlier_col]] == TRUE &
                             !is.na(df[[outlier_col]]), ]
            if (nrow(df_outlier) > 0) {
                p <- p + ggplot2::geom_point(
                    data = df_outlier,
                    ggplot2::aes(x = .data[[coord_cols[1]]],
                                 y = .data[[coord_cols[2]]]),
                    color = outlier_color, shape = 1,
                    size = point_size + 1, stroke = outlier_stroke)
            }
        }
    }

    p
}


#' Plot Per-Cluster QC Distributions
#'
#' Creates violin plots of a QC metric faceted by cluster, showing the
#' distribution within each cluster and highlighting outlier cells.
#'
#' @param spe A \linkS4class{SpatialExperiment} object with cluster labels
#'   and QC results.
#' @param cluster_col Name of the \code{colData} column with cluster labels
#'   (default "cell_cluster").
#' @param metric Name of the QC metric to plot (default "sum").
#' @param outlier_col Name of the outlier flag column (default
#'   \code{paste0(metric, "_outlier")}).
#' @param ncol Number of columns for faceting (default NULL, auto).
#'
#' @return A ggplot object.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' library(SpatialExperiment)
#' library(scuttle)
#'
#' counts <- matrix(rpois(2000, 5), nrow = 200, ncol = 10)
#' rownames(counts) <- paste0("Gene", 1:200)
#' coords <- matrix(rnorm(20), ncol = 2)
#' colnames(coords) <- c("x", "y")
#' spe <- SpatialExperiment(assays = list(counts = counts),
#'     spatialCoords = coords)
#' spe <- addPerCellQCMetrics(spe)
#' spe$cell_cluster <- rep(c("A", "B"), each = 5)
#'
#' plotClusterQC(spe, metric = "sum")
plotClusterQC <- function(spe,
                          cluster_col = "cell_cluster",
                          metric = "sum",
                          outlier_col = NULL,
                          ncol = NULL) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting. ",
             "Install with: install.packages('ggplot2')")
    }

    if (!is(spe, "SpatialExperiment")) {
        stop("'spe' must be a SpatialExperiment object.")
    }

    if (is.null(outlier_col)) {
        outlier_col <- paste0(metric, "_outlier")
    }

    cd <- as.data.frame(colData(spe))

    if (!metric %in% colnames(cd)) {
        stop("Metric '", metric, "' not found in colData.")
    }
    if (!cluster_col %in% colnames(cd)) {
        stop("Column '", cluster_col, "' not found in colData.")
    }

    has_outliers <- outlier_col %in% colnames(cd)

    p <- ggplot2::ggplot(cd, ggplot2::aes(x = .data[[cluster_col]],
                                           y = .data[[metric]])) +
        ggplot2::geom_violin(fill = "lightblue", alpha = 0.6) +
        ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste0("Per-cluster: ", metric),
                      x = "Cluster", y = metric) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                            hjust = 1))

    if (has_outliers) {
        cd_outlier <- cd[cd[[outlier_col]] == TRUE &
                         !is.na(cd[[outlier_col]]), ]
        if (nrow(cd_outlier) > 0) {
            p <- p + ggplot2::geom_jitter(
                data = cd_outlier,
                ggplot2::aes(x = .data[[cluster_col]],
                             y = .data[[metric]]),
                color = "red", alpha = 0.7, size = 0.8, width = 0.1)
        }
    }

    p
}


#' Plot Spatial Dispersion
#'
#' Creates a spatial plot colored by cluster, with artifact clusters shown in
#' grey. Renders cell boundaries as filled polygons when available, otherwise
#' falls back to points.
#'
#' @param spe A \linkS4class{SpatialExperiment} object with cluster labels
#'   and artifact flags.
#' @param cluster_col Name of the cluster label column (default
#'   "cell_cluster").
#' @param highlight_artifacts Logical. Grey out artifact clusters
#'   (default TRUE).
#' @param sample Character. Which sample to plot (default NULL, first sample).
#' @param samples Column name for sample IDs (default "sample_id").
#' @param point_size Point size for point-based rendering (default 0.5).
#' @param use_polygons Controls polygon rendering. \code{"auto"} (default)
#'   uses polygons if available, otherwise falls back to points. \code{TRUE}
#'   requires polygons. \code{FALSE} always uses points.
#' @param polygon_col Name of the \code{colData} column containing \pkg{sf}
#'   polygon geometry (default \code{"polygons"}).
#' @param polygon_border_color Border color for polygons (default
#'   \code{"grey30"}).
#' @param polygon_border_width Border line width for polygons (default 0.05).
#'
#' @return A ggplot object.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' library(SpatialExperiment)
#' library(scuttle)
#'
#' counts <- matrix(rpois(2000, 5), nrow = 200, ncol = 10)
#' rownames(counts) <- paste0("Gene", 1:200)
#' coords <- matrix(rnorm(20), ncol = 2)
#' colnames(coords) <- c("x", "y")
#' spe <- SpatialExperiment(assays = list(counts = counts),
#'     spatialCoords = coords)
#' spe <- addPerCellQCMetrics(spe)
#' spe$cell_cluster <- rep(c("A", "B"), each = 5)
#'
#' plotSpatialDispersion(spe)
plotSpatialDispersion <- function(spe,
                                  cluster_col = "cell_cluster",
                                  highlight_artifacts = TRUE,
                                  sample = NULL,
                                  samples = "sample_id",
                                  point_size = 0.5,
                                  use_polygons = "auto",
                                  polygon_col = "polygons",
                                  polygon_border_color = "grey30",
                                  polygon_border_width = 0.05) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting. ",
             "Install with: install.packages('ggplot2')")
    }

    if (!is(spe, "SpatialExperiment")) {
        stop("'spe' must be a SpatialExperiment object.")
    }

    # --- Resolve polygon vs point mode ---
    render_polygons <- FALSE
    if (identical(use_polygons, TRUE)) {
        if (!.hasPolygons(spe, polygon_col)) {
            stop("Polygon column '", polygon_col,
                 "' not found or sf package not available.")
        }
        render_polygons <- TRUE
    } else if (identical(use_polygons, "auto")) {
        render_polygons <- .hasPolygons(spe, polygon_col)
    }

    cd <- as.data.frame(colData(spe))

    # Build display_cluster column for artifact greying
    if (highlight_artifacts && "cluster_artifact_flag" %in% colnames(cd)) {
        display_cluster <- ifelse(
            cd$cluster_artifact_flag, "artifact", cd[[cluster_col]])
    } else {
        display_cluster <- cd[[cluster_col]]
    }

    # Build artifact color scale helper
    .artifact_scale <- function(display_vals) {
        if (highlight_artifacts && "artifact" %in% display_vals) {
            n_real <- length(unique(display_vals[
                display_vals != "artifact"]))
            custom_colors <- c("artifact" = "grey80")
            real_clusters <- setdiff(unique(display_vals), "artifact")
            hue_colors <- grDevices::hcl(
                h = seq(15, 375, length.out = n_real + 1)[
                    seq_len(n_real)],
                c = 100, l = 65)
            list(values = c(custom_colors,
                            stats::setNames(hue_colors, real_clusters)),
                 na.value = "grey50")
        } else {
            NULL
        }
    }

    if (render_polygons) {
        # --- Polygon-based rendering ---
        needed <- c(cluster_col, samples)
        if ("cluster_artifact_flag" %in% colnames(cd)) {
            needed <- c(needed, "cluster_artifact_flag")
        }
        sf_df <- .extractPolygonData(spe, polygon_col, needed)
        if (is.null(sf_df)) {
            render_polygons <- FALSE
        }
    }

    if (render_polygons) {
        sf_df$display_cluster <- display_cluster

        # Sample subsetting
        if (samples %in% colnames(sf_df)) {
            if (is.null(sample)) sample <- sf_df[[samples]][1]
            sf_df <- sf_df[sf_df[[samples]] == sample, ]
        }

        p <- ggplot2::ggplot(sf_df) +
            ggplot2::geom_sf(
                ggplot2::aes(fill = display_cluster),
                color = polygon_border_color,
                linewidth = polygon_border_width) +
            ggplot2::coord_sf() +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = "Spatial Dispersion", fill = "Cluster")

        scale_info <- .artifact_scale(sf_df$display_cluster)
        if (!is.null(scale_info)) {
            p <- p + ggplot2::scale_fill_manual(
                values = scale_info$values,
                na.value = scale_info$na.value)
        }
    } else {
        # --- Point-based rendering (original) ---
        coords <- as.data.frame(spatialCoords(spe))
        df <- cbind(coords, cd)
        df$display_cluster <- display_cluster
        coord_cols <- colnames(coords)

        if (samples %in% colnames(df)) {
            if (is.null(sample)) sample <- df[[samples]][1]
            df <- df[df[[samples]] == sample, ]
        }

        p <- ggplot2::ggplot(df, ggplot2::aes(
                x = .data[[coord_cols[1]]],
                y = .data[[coord_cols[2]]],
                color = display_cluster)) +
            ggplot2::geom_point(size = point_size) +
            ggplot2::coord_fixed() +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = "Spatial Dispersion", color = "Cluster")

        scale_info <- .artifact_scale(df$display_cluster)
        if (!is.null(scale_info)) {
            p <- p + ggplot2::scale_color_manual(
                values = scale_info$values,
                na.value = scale_info$na.value)
        }
    }

    p
}


#' Plot Cluster QC Summary
#'
#' Creates a multi-panel summary of cluster-level QC metrics including
#' Mahalanobis distances and spatial homogeneity scores.
#'
#' @param spe A \linkS4class{SpatialExperiment} object with cluster-level QC
#'   results from \code{\link{flagArtifactClusters}}.
#' @param cluster_col Name of the cluster label column (default
#'   "cell_cluster").
#' @param metrics QC metrics to summarize (default
#'   \code{c("sum", "detected", "subsets_mito_percent")}).
#'
#' @return A list of ggplot objects: \code{$mahalanobis} (bar plot of
#'   Mahalanobis distances), \code{$homogeneity} (bar plot of spatial
#'   homogeneity), \code{$distributions} (violin plots of selected metrics).
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' library(SpatialExperiment)
#' library(scuttle)
#'
#' counts <- matrix(rpois(2000, 5), nrow = 200, ncol = 10)
#' rownames(counts) <- paste0("Gene", 1:200)
#' coords <- matrix(rnorm(20), ncol = 2)
#' colnames(coords) <- c("x", "y")
#' spe <- SpatialExperiment(assays = list(counts = counts),
#'     spatialCoords = coords)
#' spe <- addPerCellQCMetrics(spe)
#' spe$cell_cluster <- rep(c("A", "B"), each = 5)
#'
#' plotClusterSummary(spe, metrics = "sum")
plotClusterSummary <- function(spe,
                               cluster_col = "cell_cluster",
                               metrics = c("sum", "detected",
                                           "subsets_mito_percent")) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting. ",
             "Install with: install.packages('ggplot2')")
    }

    if (!is(spe, "SpatialExperiment")) {
        stop("'spe' must be a SpatialExperiment object.")
    }

    cd <- as.data.frame(colData(spe))

    if (!cluster_col %in% colnames(cd)) {
        stop("Column '", cluster_col, "' not found in colData.")
    }

    plots <- list()

    # --- Mahalanobis distance bar plot ---
    if ("cluster_qc_mahal" %in% colnames(cd)) {
        mahal_df <- unique(cd[, c(cluster_col, "cluster_qc_mahal",
                                   "cluster_artifact_flag")])

        plots$mahalanobis <- ggplot2::ggplot(
            mahal_df,
            ggplot2::aes(x = stats::reorder(.data[[cluster_col]],
                                             cluster_qc_mahal),
                         y = cluster_qc_mahal,
                         fill = cluster_artifact_flag)) +
            ggplot2::geom_col() +
            ggplot2::scale_fill_manual(values = c("FALSE" = "steelblue",
                                                   "TRUE" = "firebrick")) +
            ggplot2::coord_flip() +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = "Cluster Mahalanobis Distance",
                          x = "Cluster", y = "Mahalanobis Distance",
                          fill = "Artifact")
    }

    # --- Spatial homogeneity bar plot ---
    if ("cluster_spatial_homogeneity" %in% colnames(cd)) {
        homo_df <- unique(cd[, c(cluster_col, "cluster_spatial_homogeneity",
                                  "cluster_artifact_flag")])

        plots$homogeneity <- ggplot2::ggplot(
            homo_df,
            ggplot2::aes(x = stats::reorder(.data[[cluster_col]],
                                             cluster_spatial_homogeneity),
                         y = cluster_spatial_homogeneity,
                         fill = cluster_artifact_flag)) +
            ggplot2::geom_col() +
            ggplot2::scale_fill_manual(values = c("FALSE" = "steelblue",
                                                   "TRUE" = "firebrick")) +
            ggplot2::coord_flip() +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = "Cluster Spatial Homogeneity",
                          x = "Cluster", y = "Neighborhood Homogeneity",
                          fill = "Artifact")
    }

    # --- Metric distributions per cluster ---
    available_metrics <- intersect(metrics, colnames(cd))
    if (length(available_metrics) > 0) {
        plots$distributions <- lapply(available_metrics, function(m) {
            plotClusterQC(spe, cluster_col = cluster_col, metric = m)
        })
        names(plots$distributions) <- available_metrics
    }

    plots
}
