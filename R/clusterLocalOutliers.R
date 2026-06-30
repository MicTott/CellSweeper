#' Within-Cluster Spatial Outlier Detection
#'
#' Detects local outlier cells within each cluster using spatially-aware
#' modified z-scores. This is Level 3 of the CellSweeper QC framework and
#' is the direct generalization of SpotSweeper to single-cell resolution:
#' spatial neighborhoods are now biologically homogeneous because they are
#' restricted to cells of the same cluster.
#'
#' For each cell, k-nearest neighbors in physical tissue coordinates are
#' found among cells of the same cluster. The cell's QC metric is then
#' compared to its neighbors using a modified z-score. Cells with extreme
#' z-scores are flagged as outliers.
#'
#' @param spe A \linkS4class{SpatialExperiment} object with cluster labels in
#'   \code{colData}.
#' @param cluster_col Name of the \code{colData} column containing cluster
#'   labels (default "cell_cluster").
#' @param exclude_artifact_clusters Logical. If TRUE (default), skip clusters
#'   flagged as artifacts by \code{\link{flagArtifactClusters}} (i.e., cells
#'   with \code{cluster_artifact_flag == TRUE}).
#' @param metrics Character vector of transcriptomic \code{colData} columns to
#'   evaluate (default: \code{c("sum", "detected", "subsets_mito_percent")}).
#' @param morpho_metrics Character vector of morphological \code{colData}
#'   columns to evaluate. Columns not found are silently skipped (default:
#'   \code{c("Area_um", "log2AspectRatio", "log2CountArea")}, matching
#'   SpaceTrooper output).
#' @param k Number of spatial nearest neighbors within the same cluster
#'   (default 50).
#' @param min_neighbors Minimum number of same-cluster spatial neighbors
#'   required for spatial QC (default 10). Cells with fewer neighbors fall
#'   back to global within-cluster z-scores.
#' @param z_threshold Z-score threshold for outlier detection (default 3).
#' @param direction Direction of outlier detection. Either a single string
#'   applied to all metrics (\code{"lower"}, \code{"higher"}, or
#'   \code{"both"}), or a named character vector mapping metric names to
#'   directions. Default: \code{"both"}.
#' @param log_transform Logical. If TRUE (default), log1p-transform skewed
#'   metrics before computing z-scores.
#' @param combine Method for combining outlier flags across metrics.
#'   \code{"union"} (default) flags a cell if outlier in ANY metric.
#'   \code{"intersection"} flags only if outlier in ALL metrics.
#' @param robust_z Scaling convention for the modified z-score. \code{"mad"}
#'   (default) uses the standard scaled MAD (\code{mad(x)} with
#'   \code{constant = 1.4826}) as a robust SD estimate, i.e.,
#'   \code{z = (x - median) / mad(x)}; threshold \code{z_threshold = 3} then
#'   corresponds to roughly 3 SD-equivalents. \code{"iglewicz"} uses the
#'   Iglewicz-Hoaglin modified z-score
#'   \code{z = 0.6745 * (x - median) / mad(x, constant = 1)};
#'   a common threshold is 3.5.
#' @param samples Column name in \code{colData} for sample IDs
#'   (default "sample_id").
#' @param BPPARAM A \code{BiocParallelParam} object for parallel
#'   processing across clusters (default NULL = serial). Set to e.g.
#'   \code{BiocParallel::MulticoreParam(4)} for 4-core parallelism.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the
#'   nearest-neighbor algorithm (default NULL = exact). Set to
#'   \code{BiocNeighbors::HnswParam()} for approximate NN on large datasets.
#'
#' @return A \linkS4class{SpatialExperiment} object with the following columns
#'   added to \code{colData} for each metric: \code{{metric}_zscore} (local
#'   z-score) and \code{{metric}_outlier} (logical flag). Also adds
#'   \code{cellsweeper_outlier} (combined flag across all metrics).
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom BiocNeighbors findKNN
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom stats median mad
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' library(CellSweeper)
#' library(SpatialExperiment)
#' library(scuttle)
#'
#' # Create a small example with cluster labels
#' counts <- matrix(rpois(2000, lambda = 5), nrow = 200, ncol = 10)
#' rownames(counts) <- paste0("Gene", 1:200)
#' colnames(counts) <- paste0("cell_", 1:10)
#' coords <- matrix(rnorm(20), ncol = 2)
#' colnames(coords) <- c("x", "y")
#'
#' spe <- SpatialExperiment(
#'     assays = list(counts = counts),
#'     spatialCoords = coords
#' )
#' spe <- addPerCellQCMetrics(spe)
#' spe$cell_cluster <- rep("A", 10)
#'
#' spe <- clusterLocalOutliers(spe, metrics = c("sum"),
#'     morpho_metrics = NULL, k = 5)
clusterLocalOutliers <- function(spe,
                                 cluster_col = "cell_cluster",
                                 exclude_artifact_clusters = TRUE,
                                 metrics = c("sum", "detected",
                                     "subsets_mito_percent"),
                                 morpho_metrics = c("Area_um",
                                     "log2AspectRatio",
                                     "log2CountArea"),
                                 k = 50,
                                 min_neighbors = 10,
                                 z_threshold = 3,
                                 direction = "both",
                                 log_transform = TRUE,
                                 combine = c("union", "intersection"),
                                 robust_z = c("mad", "iglewicz"),
                                 samples = "sample_id",
                                 BPPARAM = NULL,
                                 BNPARAM = NULL) {

    combine  <- match.arg(combine)
    robust_z <- match.arg(robust_z)

    # --- Input validation ---
    if (!is(spe, "SpatialExperiment")) {
        stop("'spe' must be a SpatialExperiment object.")
    }
    if (!cluster_col %in% colnames(colData(spe))) {
        stop("Column '", cluster_col, "' not found in colData. ",
            "Run clusterCellTypes() first.")
    }

    cd <- colData(spe)
    coords <- spatialCoords(spe)
    labels <- as.character(cd[[cluster_col]])
    n_cells <- ncol(spe)

    # Filter metrics to those present in colData
    available_metrics <- intersect(metrics, colnames(cd))
    missing_metrics <- setdiff(metrics, colnames(cd))
    if (length(missing_metrics) > 0) {
        message("clusterLocalOutliers: transcriptomic metrics not found: ",
            paste(missing_metrics, collapse = ", "))
    }

    # Filter morpho_metrics (these are expected to often be missing)
    if (!is.null(morpho_metrics)) {
        available_morpho <- intersect(morpho_metrics, colnames(cd))
        missing_morpho <- setdiff(morpho_metrics, colnames(cd))
        if (length(missing_morpho) > 0) {
            message("clusterLocalOutliers: morpho metrics not found: ",
                paste(missing_morpho, collapse = ", "))
        }
    } else {
        available_morpho <- character(0)
    }

    all_metrics <- c(available_metrics, available_morpho)
    if (length(all_metrics) == 0) {
        stop("No valid metrics found in colData.")
    }

    # --- Resolve direction per metric ---
    # Default smart directions
    default_directions <- c(
        sum = "lower", detected = "lower",
        subsets_mito_percent = "higher",
        Area_um = "both", log2AspectRatio = "both",
        log2CountArea = "lower",
        altexps_NegPrb_percent = "higher"
    )

    if (length(direction) == 1 && !is.null(names(direction))) {
        # Single named element <U+2014> treat as single direction for all
        metric_directions <- rep(direction, length(all_metrics))
        names(metric_directions) <- all_metrics
    } else if (length(direction) == 1) {
        metric_directions <- rep(direction, length(all_metrics))
        names(metric_directions) <- all_metrics
    } else {
        # Named vector: use provided, fall back to defaults
        metric_directions <- vapply(all_metrics, function(m) {
            if (m %in% names(direction)) direction[m]
            else if (m %in% names(default_directions)) default_directions[m]
            else "both"
        }, character(1))
    }

    # Metrics that should be log-transformed (right-skewed)
    # Note: log2AspectRatio, log2CountArea are already log-scaled
    log_candidates <- c("sum", "detected", "Area_um")

    # --- Initialize output columns ---
    for (m in all_metrics) {
        cd[[paste0(m, "_zscore")]] <- NA_real_
        cd[[paste0(m, "_outlier")]] <- FALSE
    }

    # --- Handle sample_id ---
    if (samples %in% colnames(cd)) {
        sample_vec <- as.character(cd[[samples]])
    } else {
        sample_vec <- rep("all", n_cells)
    }
    sample_ids <- unique(sample_vec)

    # --- Per-task worker (operates on pre-sliced data) ---
    # Each `task` carries its own `coords`, `vals`, and `indices`; the
    # full coord matrix and metric vectors do NOT travel into worker
    # processes under bplapply.
    .process_cluster <- function(task, all_metrics, metric_directions,
                                 log_transform, log_candidates,
                                 k, min_neighbors, z_threshold,
                                 robust_z, BNPARAM) {
        cluster_indices <- task$indices
        n_cluster <- length(cluster_indices)
        if (n_cluster < 3) return(NULL)

        # Compute kNN within this cluster
        k_actual <- min(k, n_cluster - 1)
        fknn_args <- list(X = task$coords, k = k_actual,
            warn.ties = FALSE)
        if (!is.null(BNPARAM)) fknn_args$BNPARAM <- BNPARAM
        nn <- do.call(BiocNeighbors::findKNN, fknn_args)$index

        result <- list()
        for (m in all_metrics) {
            values <- task$vals[[m]]

            if (log_transform && m %in% log_candidates) {
                values <- log1p(pmax(values, 0))
            }

            dir_m <- metric_directions[m]

            z_scores <- vapply(seq_len(n_cluster), function(i) {
                neighbor_idx <- nn[i, ]
                neighbor_idx <- neighbor_idx[neighbor_idx != 0]
                if (length(neighbor_idx) < min_neighbors) {
                    neighborhood <- values[-i]
                } else {
                    neighborhood <- values[neighbor_idx]
                }
                if (length(neighborhood) == 0) return(0)
                med <- stats::median(neighborhood)
                if (robust_z == "iglewicz") {
                    # Iglewicz-Hoaglin: MAD_raw (constant = 1), prefactor 0.6745
                    mad_val <- stats::mad(neighborhood, center = med,
                        constant = 1)
                    if (mad_val == 0) return(0)
                    z <- 0.6745 * (values[i] - med) / mad_val
                } else {
                    # Robust z using scaled MAD (constant = 1.4826) as
                    # an estimate of SD; threshold ~ Z-score-equivalent
                    mad_val <- stats::mad(neighborhood, center = med)
                    if (mad_val == 0) return(0)
                    z <- (values[i] - med) / mad_val
                }
                if (!is.finite(z)) 0 else z
            }, numeric(1))

            outlier_flags <- switch(dir_m,
                lower  = z_scores < -z_threshold,
                higher = z_scores >  z_threshold,
                both   = abs(z_scores) > z_threshold,
                abs(z_scores) > z_threshold
            )

            result[[paste0(m, "_zscore")]]  <- z_scores
            result[[paste0(m, "_outlier")]] <- outlier_flags
        }
        list(indices = cluster_indices, results = result)
    }

    # Build task list <U+2014> pre-slice coords and metric vectors per task so
    # bplapply does not serialize the entire dataset N times.
    cd_vals <- lapply(all_metrics, function(m) as.numeric(cd[[m]]))
    names(cd_vals) <- all_metrics

    tasks <- list()
    for (sid in sample_ids) {
        sample_mask    <- sample_vec == sid
        sample_labels  <- labels[sample_mask]
        sample_indices <- which(sample_mask)

        for (cl in unique(sample_labels)) {
            cluster_indices <- sample_indices[sample_labels == cl]
            if (exclude_artifact_clusters &&
                "cluster_artifact_flag" %in% colnames(cd)) {
                if (length(cluster_indices) > 0 &&
                    all(cd$cluster_artifact_flag[cluster_indices])) {
                    next
                }
            }
            tasks <- c(tasks, list(list(
                indices = cluster_indices,
                coords  = coords[cluster_indices, , drop = FALSE],
                vals    = lapply(cd_vals,
                    function(v) v[cluster_indices])
            )))
        }
    }

    # Run tasks (parallel or serial)
    .run_args <- list(
        all_metrics = all_metrics,
        metric_directions = metric_directions,
        log_transform = log_transform,
        log_candidates = log_candidates,
        k = k, min_neighbors = min_neighbors,
        z_threshold = z_threshold,
        robust_z = robust_z,
        BNPARAM = BNPARAM
    )
    if (!is.null(BPPARAM)) {
        results <- do.call(BiocParallel::bplapply,
            c(list(X = tasks, FUN = .process_cluster), .run_args,
                list(BPPARAM = BPPARAM)))
    } else {
        results <- do.call(lapply,
            c(list(X = tasks, FUN = .process_cluster), .run_args))
    }

    # Collect results back into colData
    for (res in results) {
        if (is.null(res)) next
        for (col_name in names(res$results)) {
            cd[[col_name]][res$indices] <- res$results[[col_name]]
        }
    }

    # --- Combined outlier flag ---
    outlier_cols <- paste0(all_metrics, "_outlier")
    outlier_mat <- as.matrix(
        as.data.frame(cd[, outlier_cols, drop = FALSE]))

    # Handle NAs (cells in artifact clusters or tiny clusters)
    outlier_mat[is.na(outlier_mat)] <- FALSE

    if (combine == "union") {
        cd$cellsweeper_outlier <- apply(outlier_mat, 1, any)
    } else {
        cd$cellsweeper_outlier <- apply(outlier_mat, 1, all)
    }

    colData(spe) <- cd

    n_outliers <- sum(cd$cellsweeper_outlier, na.rm = TRUE)
    message("clusterLocalOutliers: ", n_outliers, "/", n_cells,
        " cells flagged as outliers (",
        round(100 * n_outliers / n_cells, 1), "%)")

    spe
}
