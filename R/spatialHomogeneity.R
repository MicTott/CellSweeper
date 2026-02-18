#' Compute Neighborhood Homogeneity
#'
#' For each cell, computes the fraction of its k spatial nearest neighbors
#' that share its cluster label. This measures how spatially coherent a
#' cluster is — real cell types tend to be spatially organized, while
#' artifact clusters (driven by shared damage signatures) tend to be
#' spatially scattered.
#'
#' @param spe A \linkS4class{SpatialExperiment} object with cluster labels
#'   in \code{colData}.
#' @param cluster_col Name of the \code{colData} column containing cluster
#'   labels (default "cell_cluster").
#' @param k Number of spatial nearest neighbors to use (default 20).
#' @param samples Column name in \code{colData} for sample IDs
#'   (default "sample_id"). Spatial neighbors are computed within samples.
#'
#' @return A \linkS4class{SpatialExperiment} object with
#'   \code{neighborhood_homogeneity} added to \code{colData}.
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom BiocNeighbors findKNN
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
#' spe$cell_cluster <- rep(c("A", "B"), each = 5)
#'
#' spe <- computeNeighborhoodHomogeneity(spe, k = 3)
computeNeighborhoodHomogeneity <- function(spe,
                                           cluster_col = "cell_cluster",
                                           k = 20,
                                           samples = "sample_id") {

    if (!is(spe, "SpatialExperiment")) {
        stop("'spe' must be a SpatialExperiment object.")
    }
    if (!cluster_col %in% colnames(colData(spe))) {
        stop("Column '", cluster_col, "' not found in colData.")
    }

    cd <- colData(spe)
    coords <- spatialCoords(spe)
    labels <- as.character(cd[[cluster_col]])
    n_cells <- ncol(spe)

    homogeneity <- rep(NA_real_, n_cells)

    # Process per-sample (spatial neighbors must be within same tissue section)
    if (samples %in% colnames(cd)) {
        sample_ids <- unique(cd[[samples]])
    } else {
        cd[["_tmp_sample"]] <- "all"
        sample_ids <- "all"
        samples <- "_tmp_sample"
    }

    for (sid in sample_ids) {
        idx <- which(cd[[samples]] == sid)
        if (length(idx) < 2) next

        k_actual <- min(k, length(idx) - 1)
        nn <- BiocNeighbors::findKNN(coords[idx, , drop = FALSE],
                                     k = k_actual,
                                     warn.ties = FALSE)$index

        sample_labels <- labels[idx]

        # For each cell, fraction of neighbors sharing its label
        for (i in seq_along(idx)) {
            neighbor_idx <- nn[i, ]
            neighbor_idx <- neighbor_idx[neighbor_idx != 0]
            if (length(neighbor_idx) == 0) {
                homogeneity[idx[i]] <- NA_real_
            } else {
                homogeneity[idx[i]] <-
                    mean(sample_labels[neighbor_idx] == sample_labels[i])
            }
        }
    }

    colData(spe)$neighborhood_homogeneity <- homogeneity
    spe
}


#' Permutation Test for Spatial Coherence
#'
#' Computes a permutation null for neighborhood homogeneity to determine
#' whether each cluster's spatial coherence is greater than expected by
#' chance. Large clusters have higher baseline homogeneity simply due to
#' their size; the permutation test accounts for this.
#'
#' @param spe A \linkS4class{SpatialExperiment} object with cluster labels
#'   and \code{neighborhood_homogeneity} in \code{colData} (from
#'   \code{\link{computeNeighborhoodHomogeneity}}).
#' @param cluster_col Name of the \code{colData} column containing cluster
#'   labels (default "cell_cluster").
#' @param k Number of spatial nearest neighbors (default 20).
#' @param n_permutations Number of permutations for the null distribution
#'   (default 500).
#' @param samples Column name in \code{colData} for sample IDs
#'   (default "sample_id").
#' @param seed Random seed for reproducibility (default 42).
#'
#' @return A data.frame with one row per cluster and columns:
#'   \code{cluster}, \code{observed_homogeneity}, \code{null_mean},
#'   \code{null_sd}, \code{z_score}, \code{p_value}.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom BiocNeighbors findKNN
#' @importFrom stats sd
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' library(SpatialExperiment)
#'
#' counts <- matrix(rpois(2000, 5), nrow = 200, ncol = 10)
#' rownames(counts) <- paste0("Gene", 1:200)
#' coords <- matrix(rnorm(20), ncol = 2)
#' colnames(coords) <- c("x", "y")
#' spe <- SpatialExperiment(assays = list(counts = counts),
#'     spatialCoords = coords)
#' spe$cell_cluster <- rep(c("A", "B"), each = 5)
#' spe <- computeNeighborhoodHomogeneity(spe, k = 3)
#'
#' results <- permuteHomogeneity(spe, k = 3, n_permutations = 10)
permuteHomogeneity <- function(spe,
                               cluster_col = "cell_cluster",
                               k = 20,
                               n_permutations = 500,
                               samples = "sample_id",
                               seed = 42) {

    if (!is(spe, "SpatialExperiment")) {
        stop("'spe' must be a SpatialExperiment object.")
    }
    if (!cluster_col %in% colnames(colData(spe))) {
        stop("Column '", cluster_col, "' not found in colData.")
    }

    cd <- colData(spe)
    coords <- spatialCoords(spe)
    labels <- as.character(cd[[cluster_col]])
    unique_clusters <- unique(labels)
    n_cells <- ncol(spe)

    # Handle sample_id
    if (samples %in% colnames(cd)) {
        sample_vec <- as.character(cd[[samples]])
    } else {
        sample_vec <- rep("all", n_cells)
    }
    sample_ids <- unique(sample_vec)

    # --- Pre-compute kNN index matrices per sample ---
    nn_list <- list()
    idx_list <- list()
    for (sid in sample_ids) {
        idx <- which(sample_vec == sid)
        idx_list[[sid]] <- idx
        if (length(idx) < 2) {
            nn_list[[sid]] <- NULL
            next
        }
        k_actual <- min(k, length(idx) - 1)
        nn_list[[sid]] <- BiocNeighbors::findKNN(
            coords[idx, , drop = FALSE],
            k = k_actual, warn.ties = FALSE)$index
    }

    # --- Observed per-cluster mean homogeneity ---
    .compute_cluster_homogeneity <- function(labs) {
        # Vectorized computation across all samples
        homo <- rep(NA_real_, n_cells)
        for (sid in sample_ids) {
            idx <- idx_list[[sid]]
            nn <- nn_list[[sid]]
            if (is.null(nn)) next
            sample_labs <- labs[idx]
            # Vectorized: compare each cell's label to its neighbors' labels
            neighbor_labels <- matrix(sample_labs[nn], nrow = nrow(nn))
            homo[idx] <- rowMeans(neighbor_labels == sample_labs,
                                  na.rm = TRUE)
        }
        # Per-cluster means
        vapply(unique_clusters, function(cl) {
            mean(homo[labs == cl], na.rm = TRUE)
        }, numeric(1))
    }

    observed <- .compute_cluster_homogeneity(labels)

    # --- Permutation null ---
    null_matrix <- matrix(NA_real_, nrow = length(unique_clusters),
                          ncol = n_permutations)
    rownames(null_matrix) <- unique_clusters

    for (p in seq_len(n_permutations)) {
        # Permute cluster labels (within each sample to respect sample structure)
        perm_labels <- labels
        for (sid in sample_ids) {
            idx <- idx_list[[sid]]
            perm_labels[idx] <- sample(perm_labels[idx])
        }
        null_matrix[, p] <- .compute_cluster_homogeneity(perm_labels)
    }

    # --- Compute z-scores and p-values ---
    null_means <- rowMeans(null_matrix)
    null_sds <- apply(null_matrix, 1, sd)

    z_scores <- ifelse(null_sds == 0, 0,
                        (observed - null_means) / null_sds)

    # One-sided p-value: proportion of permutations >= observed
    p_values <- vapply(seq_along(unique_clusters), function(i) {
        (sum(null_matrix[i, ] >= observed[i]) + 1) / (n_permutations + 1)
    }, numeric(1))

    data.frame(
        cluster = unique_clusters,
        observed_homogeneity = observed,
        null_mean = null_means,
        null_sd = null_sds,
        z_score = z_scores,
        p_value = p_values,
        stringsAsFactors = FALSE
    )
}
