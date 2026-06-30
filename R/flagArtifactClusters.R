#' Flag Artifact Clusters
#'
#' Identifies entire clusters that are likely technical artifacts rather than
#' real cell types. This is Level 2b+c of the CellSweeper QC framework.
#' Uses two complementary approaches:
#'
#' \enumerate{
#'   \item \strong{Pseudobulk multivariate outlier detection:} Computes
#'     per-cluster summary statistics (medians of QC metrics) and identifies
#'     clusters with extreme Mahalanobis distance from the center.
#'   \item \strong{Spatial dispersion analysis:} Real cell types are spatially
#'     organized in tissue; artifact clusters tend to be spatially scattered.
#'     A permutation test assesses whether each cluster is more spatially
#'     coherent than expected by chance.
#' }
#'
#' A cluster that is BOTH a QC outlier AND spatially incoherent is flagged as
#' an artifact. A low-expression cluster with strong spatial coherence is
#' likely a legitimate cell type and is NOT flagged.
#'
#' @param spe A \linkS4class{SpatialExperiment} object with cluster labels in
#'   \code{colData}.
#' @param cluster_col Name of the \code{colData} column containing cluster
#'   labels (default "cell_cluster").
#' @param metrics Character vector of \code{colData} column names to use for
#'   multivariate outlier detection. Columns not found are skipped with a
#'   message. Default: \code{c("sum", "detected", "subsets_mito_percent",
#'   "Area_um", "log2CountArea")} (matching SpaceTrooper output).
#' @param mahal_threshold Threshold for Mahalanobis distance outlier detection.
#'   Use "auto" (default) for chi-squared based threshold at p = 0.05, or a
#'   numeric value.
#' @param covariance_method Method for covariance estimation. "mcd" (default)
#'   uses minimum covariance determinant via \code{robustbase::covMcd()} if
#'   available, otherwise falls back to \code{stats::cov()}.
#' @param spatial_k Number of spatial nearest neighbors for homogeneity
#'   computation (default 20).
#' @param n_permutations Number of permutations for spatial coherence null
#'   distribution (default 500).
#' @param dispersion_threshold P-value threshold for spatial incoherence
#'   (default 0.05). Clusters with p-value above this are considered spatially
#'   incoherent.
#' @param require_both Logical. If TRUE (default), a cluster must be BOTH a QC
#'   outlier AND spatially incoherent to be flagged. If FALSE, either condition
#'   suffices.
#' @param samples Column name in \code{colData} for sample IDs
#'   (default "sample_id").
#' @param BPPARAM A \code{BiocParallelParam} object for parallel
#'   permutation testing (default NULL = serial). Passed to
#'   \code{\link{permuteHomogeneity}}.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object for nearest-neighbor
#'   search (default NULL = exact). Passed to
#'   \code{\link{computeNeighborhoodHomogeneity}} and
#'   \code{\link{permuteHomogeneity}}.
#' @param seed Random seed for reproducibility (default 42).
#'
#' @return A \linkS4class{SpatialExperiment} object with the following columns
#'   added to \code{colData}: \code{cluster_qc_mahal} (Mahalanobis distance of
#'   cell's cluster), \code{cluster_spatial_homogeneity} (neighborhood
#'   homogeneity), \code{cluster_spatial_pvalue} (permutation p-value),
#'   \code{cluster_artifact_flag} (logical, TRUE = artifact cluster).
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom stats cov mahalanobis median qchisq prcomp
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
#' spe$cell_cluster <- rep(c("A", "B"), each = 5)
#'
#' # Note: with only 2 clusters, Mahalanobis is skipped
#' spe <- flagArtifactClusters(spe, metrics = c("sum", "detected"),
#'     n_permutations = 10)
flagArtifactClusters <- function(spe,
                                 cluster_col = "cell_cluster",
                                 metrics = c("sum", "detected",
                                     "subsets_mito_percent",
                                     "Area_um",
                                     "log2CountArea"),
                                 mahal_threshold = "auto",
                                 covariance_method = "mcd",
                                 spatial_k = 20,
                                 n_permutations = 500,
                                 dispersion_threshold = 0.05,
                                 require_both = TRUE,
                                 samples = "sample_id",
                                 BPPARAM = NULL,
                                 BNPARAM = NULL,
                                 seed = 42) {
    # --- Input validation ---
    if (!is(spe, "SpatialExperiment")) {
        stop("'spe' must be a SpatialExperiment object.")
    }
    if (!cluster_col %in% colnames(colData(spe))) {
        stop("Column '", cluster_col, "' not found in colData. ",
            "Run clusterCellTypes() first.")
    }

    cd <- colData(spe)
    labels <- as.character(cd[[cluster_col]])
    n_cells <- ncol(spe)

    # Handle sample_id
    if (samples %in% colnames(cd)) {
        sample_vec <- as.character(cd[[samples]])
    } else {
        sample_vec <- rep("all", n_cells)
    }
    sample_ids <- unique(sample_vec)

    # Filter metrics to those present in colData
    available_metrics <- intersect(metrics, colnames(cd))
    missing_metrics <- setdiff(metrics, colnames(cd))
    if (length(missing_metrics) > 0) {
        message("flagArtifactClusters: metrics not found, skipping: ",
            paste(missing_metrics, collapse = ", "))
    }
    if (length(available_metrics) == 0) {
        stop("No valid metrics found in colData.")
    }

    # ================================================================
    # Part 1: Per-Sample Pseudobulk Multivariate Outlier Detection
    # ================================================================

    cell_mahal         <- rep(NA_real_, n_cells)
    cell_mahal_outlier <- rep(FALSE, n_cells)
    # Track whether Mahalanobis was actually computable for each sample.
    # When it could NOT be computed (< 3 clusters in that sample), we
    # MUST NOT silently AND it against the spatial flag below: that would
    # disable artifact detection entirely. We fall back to the spatial
    # flag alone for those cells.
    cell_mahal_skipped <- rep(FALSE, n_cells)

    for (sid in sample_ids) {
        sample_idx <- which(sample_vec == sid)
        sample_labels <- labels[sample_idx]
        sample_clusters <- unique(sample_labels)
        n_cl <- length(sample_clusters)

        # Compute per-cluster medians for this sample
        cluster_summaries <- do.call(rbind, lapply(sample_clusters,
            function(cl) {
                mask <- sample_labels == cl
                vapply(available_metrics, function(m) {
                    median(as.numeric(cd[[m]][sample_idx[mask]]),
                        na.rm = TRUE)
                }, numeric(1))
            }))
        rownames(cluster_summaries) <- sample_clusters

        mahal_distances <- rep(NA_real_, n_cl)
        names(mahal_distances) <- sample_clusters
        mahal_outlier <- rep(FALSE, n_cl)
        names(mahal_outlier) <- sample_clusters

        if (n_cl < 3) {
            cell_mahal_skipped[sample_idx] <- TRUE
            if (length(sample_ids) == 1) {
                warning("Fewer than 3 clusters (", n_cl,
                    "). Skipping Mahalanobis distance-based outlier ",
                    "detection.")
            } else {
                message("flagArtifactClusters: sample '", sid, "' has ",
                    n_cl, " cluster(s), skipping Mahalanobis.")
            }
        } else {
            n_metrics <- ncol(cluster_summaries)

            if (n_metrics >= n_cl) {
                message("flagArtifactClusters: sample '", sid,
                    "': more metrics (", n_metrics,
                    ") than clusters (", n_cl,
                    "). Reducing via PCA.")
                pca_result <- stats::prcomp(cluster_summaries, scale. = TRUE)
                n_pcs <- min(n_cl - 1, n_metrics)
                cluster_summaries_use <- pca_result$x[, seq_len(n_pcs),
                    drop = FALSE]
            } else {
                cluster_summaries_use <- scale(cluster_summaries)
            }

            cov_mat <- tryCatch(
                {
                    if (covariance_method == "mcd" &&
                        requireNamespace("robustbase", quietly = TRUE)) {
                        mcd <- robustbase::covMcd(cluster_summaries_use)
                        list(center = mcd$center, cov = mcd$cov)
                    } else {
                        if (covariance_method == "mcd") {
                            message("flagArtifactClusters: robustbase not ",
                                "available, using standard covariance.")
                        }
                        list(center = colMeans(cluster_summaries_use),
                            cov = cov(cluster_summaries_use))
                    }
                },
                error = function(e) {
                    message("flagArtifactClusters: covariance estimation ",
                        "failed, using standard covariance. Reason: ",
                        e$message)
                    list(center = colMeans(cluster_summaries_use),
                        cov = cov(cluster_summaries_use))
                })

            mahal_distances <- tryCatch(
                {
                    mahalanobis(cluster_summaries_use,
                        center = cov_mat$center,
                        cov = cov_mat$cov)
                },
                error = function(e) {
                    message("flagArtifactClusters: Mahalanobis computation ",
                        "failed. Reason: ", e$message)
                    rep(0, n_cl)
                })
            names(mahal_distances) <- sample_clusters

            if (identical(mahal_threshold, "auto")) {
                df <- ncol(cluster_summaries_use)
                threshold <- qchisq(0.95, df = df)
            } else {
                threshold <- as.numeric(mahal_threshold)
            }

            mahal_outlier <- mahal_distances > threshold
        }

        # Map cluster-level Mahalanobis to per-cell for this sample
        for (cl in sample_clusters) {
            cell_mask <- sample_idx[sample_labels == cl]
            cell_mahal[cell_mask] <- mahal_distances[cl]
            cell_mahal_outlier[cell_mask] <- mahal_outlier[cl]
        }
    }

    # ================================================================
    # Part 2: Spatial Dispersion Analysis
    # ================================================================

    spe <- computeNeighborhoodHomogeneity(spe, cluster_col = cluster_col,
        k = spatial_k, samples = samples,
        BNPARAM = BNPARAM)

    perm_results <- permuteHomogeneity(spe, cluster_col = cluster_col,
        k = spatial_k,
        n_permutations = n_permutations,
        samples = samples,
        BPPARAM = BPPARAM,
        BNPARAM = BNPARAM,
        seed = seed)

    # Map spatial results to per-cell vectors
    cell_pvalue <- rep(NA_real_, n_cells)
    cell_homo <- rep(NA_real_, n_cells)
    cell_spatial_incoherent <- rep(FALSE, n_cells)

    # `p_value` from permuteHomogeneity() is the UPPER-tail probability
    # P(null homogeneity >= observed). Spatially coherent clusters have
    # observed homogeneity well above the null, hence a small p; the
    # condition `p > dispersion_threshold` therefore marks clusters that
    # FAILED to demonstrate spatial coherence and are candidate artifacts.
    for (i in seq_len(nrow(perm_results))) {
        sid  <- perm_results$sample[i]
        cl   <- perm_results$cluster[i]
        mask <- which(sample_vec == sid & labels == cl)
        cell_pvalue[mask] <- perm_results$p_value[i]
        cell_homo[mask]   <- perm_results$observed_homogeneity[i]
        cell_spatial_incoherent[mask] <-
            perm_results$p_value[i] > dispersion_threshold
    }

    # ================================================================
    # Part 3: Combine Evidence (per-cell)
    # ================================================================

    cell_artifact <- if (require_both) {
        # Where Mahalanobis is available, require BOTH conditions.
        # Where Mahalanobis was skipped (< 3 clusters in a sample),
        # fall back to the spatial flag alone <U+2014> otherwise the
        # `& FALSE` would silently disable artifact detection.
        ifelse(cell_mahal_skipped,
            cell_spatial_incoherent,
            cell_mahal_outlier & cell_spatial_incoherent)
    } else {
        cell_mahal_outlier | cell_spatial_incoherent
    }

    # ================================================================
    # Store results in colData
    # ================================================================

    colData(spe)$cluster_qc_mahal <- cell_mahal
    colData(spe)$cluster_spatial_homogeneity <- cell_homo
    colData(spe)$cluster_spatial_pvalue <- cell_pvalue
    colData(spe)$cluster_artifact_flag <- cell_artifact

    n_flagged_clusters <- sum(cell_artifact)
    message("flagArtifactClusters: ", n_flagged_clusters, "/", n_cells,
        " cell(s) flagged as artifacts")

    spe
}
