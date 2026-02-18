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
#'                             n_permutations = 10)
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
    unique_clusters <- unique(labels)
    n_clusters <- length(unique_clusters)

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
    # Part 1: Pseudobulk Multivariate Outlier Detection
    # ================================================================

    # Compute per-cluster medians for each metric
    cluster_summaries <- do.call(rbind, lapply(unique_clusters, function(cl) {
        mask <- labels == cl
        vapply(available_metrics, function(m) {
            median(as.numeric(cd[[m]][mask]), na.rm = TRUE)
        }, numeric(1))
    }))
    rownames(cluster_summaries) <- unique_clusters

    # Mahalanobis distance
    mahal_distances <- rep(NA_real_, n_clusters)
    names(mahal_distances) <- unique_clusters
    mahal_outlier <- rep(FALSE, n_clusters)
    names(mahal_outlier) <- unique_clusters

    if (n_clusters < 3) {
        warning("Fewer than 3 clusters (", n_clusters,
                "). Skipping Mahalanobis distance-based outlier detection.")
    } else {
        n_metrics <- ncol(cluster_summaries)

        # If more metrics than clusters, reduce dimensionality
        if (n_metrics >= n_clusters) {
            message("flagArtifactClusters: more metrics (", n_metrics,
                    ") than clusters (", n_clusters,
                    "). Reducing via PCA.")
            pca_result <- stats::prcomp(cluster_summaries, scale. = TRUE)
            n_pcs <- min(n_clusters - 1, n_metrics)
            cluster_summaries_use <- pca_result$x[, seq_len(n_pcs),
                                                    drop = FALSE]
        } else {
            cluster_summaries_use <- scale(cluster_summaries)
        }

        # Covariance estimation
        cov_mat <- tryCatch({
            if (covariance_method == "mcd" &&
                requireNamespace("robustbase", quietly = TRUE)) {
                mcd <- robustbase::covMcd(cluster_summaries_use)
                list(center = mcd$center, cov = mcd$cov)
            } else {
                if (covariance_method == "mcd") {
                    message("flagArtifactClusters: robustbase not available, ",
                            "using standard covariance.")
                }
                list(center = colMeans(cluster_summaries_use),
                     cov = cov(cluster_summaries_use))
            }
        }, error = function(e) {
            message("flagArtifactClusters: covariance estimation failed, ",
                    "using standard covariance. Reason: ", e$message)
            list(center = colMeans(cluster_summaries_use),
                 cov = cov(cluster_summaries_use))
        })

        # Compute Mahalanobis distances
        mahal_distances <- tryCatch({
            mahalanobis(cluster_summaries_use,
                        center = cov_mat$center,
                        cov = cov_mat$cov)
        }, error = function(e) {
            message("flagArtifactClusters: Mahalanobis computation failed. ",
                    "Reason: ", e$message)
            rep(0, n_clusters)
        })
        names(mahal_distances) <- unique_clusters

        # Threshold
        if (identical(mahal_threshold, "auto")) {
            df <- ncol(cluster_summaries_use)
            threshold <- qchisq(0.95, df = df)
        } else {
            threshold <- as.numeric(mahal_threshold)
        }

        mahal_outlier <- mahal_distances > threshold
    }

    # ================================================================
    # Part 2: Spatial Dispersion Analysis
    # ================================================================

    spe <- computeNeighborhoodHomogeneity(spe, cluster_col = cluster_col,
                                          k = spatial_k, samples = samples)

    perm_results <- permuteHomogeneity(spe, cluster_col = cluster_col,
                                       k = spatial_k,
                                       n_permutations = n_permutations,
                                       samples = samples, seed = seed)

    # Spatially incoherent: high p-value means NOT more coherent than random
    spatial_incoherent <- rep(FALSE, n_clusters)
    names(spatial_incoherent) <- unique_clusters
    for (i in seq_len(nrow(perm_results))) {
        cl <- perm_results$cluster[i]
        spatial_incoherent[cl] <- perm_results$p_value[i] > dispersion_threshold
    }

    # ================================================================
    # Part 3: Combine Evidence
    # ================================================================

    artifact_flag <- if (require_both) {
        mahal_outlier & spatial_incoherent
    } else {
        mahal_outlier | spatial_incoherent
    }

    # ================================================================
    # Store results in colData (per-cell, inherited from cluster)
    # ================================================================

    # Map cluster-level results to per-cell
    cluster_to_mahal <- mahal_distances[labels]
    cluster_to_pvalue <- rep(NA_real_, ncol(spe))
    cluster_to_homo <- rep(NA_real_, ncol(spe))
    cluster_to_artifact <- rep(FALSE, ncol(spe))

    for (i in seq_len(nrow(perm_results))) {
        cl <- perm_results$cluster[i]
        mask <- labels == cl
        cluster_to_pvalue[mask] <- perm_results$p_value[i]
        cluster_to_homo[mask] <- perm_results$observed_homogeneity[i]
        cluster_to_artifact[mask] <- artifact_flag[cl]
    }

    colData(spe)$cluster_qc_mahal <- as.numeric(cluster_to_mahal)
    colData(spe)$cluster_spatial_homogeneity <- cluster_to_homo
    colData(spe)$cluster_spatial_pvalue <- cluster_to_pvalue
    colData(spe)$cluster_artifact_flag <- cluster_to_artifact

    n_flagged <- sum(artifact_flag)
    message("flagArtifactClusters: ", n_flagged, "/", n_clusters,
            " cluster(s) flagged as artifacts")

    spe
}
