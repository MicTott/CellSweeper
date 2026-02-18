#' Cluster Cells into Coarse Cell Types
#'
#' Performs coarse clustering of cells in gene expression space to group
#' biologically similar cells. This is Level 2a of the CellSweeper QC
#' framework. Clustering restores the local homogeneity assumption needed
#' for within-cluster spatial outlier detection.
#'
#' The function log-normalizes expression, selects highly variable genes,
#' runs PCA, and clusters using a shared nearest-neighbor graph. Optionally,
#' QC-related covariates (library size, cell area, etc.) can be regressed
#' out of the expression matrix before PCA to prevent QC metrics from
#' driving the clustering.
#'
#' @param spe A \linkS4class{SpatialExperiment} object with a \code{counts}
#'   assay.
#' @param resolution Leiden clustering resolution (default 0.5). Lower values
#'   produce fewer, coarser clusters.
#' @param regress_covariates Logical. If TRUE (default), regress QC covariates
#'   from the log-normalized expression before PCA.
#' @param covariates Character vector of \code{colData} column names to regress
#'   out (default: \code{c("sum", "detected", "subsets_mito_percent")}).
#'   Columns not found in \code{colData} are silently skipped.
#' @param exclude_mito Logical. If TRUE (default), exclude mitochondrial genes
#'   from feature selection.
#' @param exclude_ribo Logical. If TRUE (default), exclude ribosomal genes
#'   from feature selection.
#' @param mito_pattern Regex pattern for mitochondrial genes
#'   (default "^MT-|^mt-").
#' @param ribo_pattern Regex pattern for ribosomal genes
#'   (default \code{"^RP[SL]|^Rp[sl]"}).
#' @param n_hvgs Number of highly variable genes to select (default 2000).
#' @param n_pcs Number of principal components to compute (default 30).
#' @param min_cluster_size Minimum cells per cluster (default 50). Clusters
#'   smaller than this are labeled "unassigned" and fall back to global QC
#'   in Level 3.
#' @param k.neighbors Number of nearest neighbors for the shared nearest-
#'   neighbor graph used in Leiden clustering (default 25).
#' @param cluster_col Name of the \code{colData} column to store cluster
#'   labels (default "cell_cluster").
#' @param use_existing Name of an existing \code{colData} column containing
#'   pre-computed cluster labels. If provided, clustering is skipped entirely
#'   and these labels are used. Set to NULL (default) to perform clustering.
#' @param seed Random seed for reproducibility (default 42).
#'
#' @return A \linkS4class{SpatialExperiment} object with cluster labels in
#'   \code{colData(spe)[[cluster_col]]}. PCA results are stored in
#'   \code{reducedDim(spe, "PCA")}.
#'
#' @importFrom SummarizedExperiment colData colData<- assay assay<-
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom scran clusterCells modelGeneVar getTopHVGs
#' @importFrom scater logNormCounts runPCA
#' @importFrom bluster NNGraphParam
#' @importFrom stats lm.fit model.matrix
#' @importFrom methods is
#' @importFrom Matrix t
#'
#' @export
#'
#' @examples
#' library(CellSweeper)
#' library(SpatialExperiment)
#' library(scuttle)
#'
#' # Create a small example
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
#'
#' # Use pre-existing labels instead of clustering
#' spe$my_labels <- rep(c("A", "B"), each = 5)
#' spe <- clusterCellTypes(spe, use_existing = "my_labels")
clusterCellTypes <- function(spe,
                             resolution = 0.5,
                             regress_covariates = TRUE,
                             covariates = c("sum", "detected",
                                            "subsets_mito_percent"),
                             exclude_mito = TRUE,
                             exclude_ribo = TRUE,
                             mito_pattern = "^MT-|^mt-",
                             ribo_pattern = "^RP[SL]|^Rp[sl]",
                             n_hvgs = 2000,
                             n_pcs = 30,
                             min_cluster_size = 50,
                             k.neighbors=25,
                             cluster_col = "cell_cluster",
                             use_existing = NULL,
                             seed = 42) {

    # --- Input validation ---
    if (!is(spe, "SpatialExperiment")) {
        stop("'spe' must be a SpatialExperiment object.")
    }

    if (!is.null(use_existing)) {
        if (!use_existing %in% colnames(colData(spe))) {
            stop("Column '", use_existing, "' not found in colData.")
        }
        colData(spe)[[cluster_col]] <- as.character(
            colData(spe)[[use_existing]])
        message("clusterCellTypes: using existing labels from '",
                use_existing, "'")
        return(spe)
    }

    if (ncol(spe) < min_cluster_size) {
        warning("Fewer cells (", ncol(spe), ") than min_cluster_size (",
                min_cluster_size, "). Assigning all cells to one cluster.")
        colData(spe)[[cluster_col]] <- "cluster_1"
        return(spe)
    }

    # --- Gene filtering ---
    genes_to_use <- rownames(spe)

    if (exclude_mito) {
        mito_genes <- grep(mito_pattern, genes_to_use, value = TRUE)
        if (length(mito_genes) > 0) {
            genes_to_use <- setdiff(genes_to_use, mito_genes)
            message("clusterCellTypes: excluded ", length(mito_genes),
                    " mitochondrial genes")
        }
    }

    if (exclude_ribo) {
        ribo_genes <- grep(ribo_pattern, genes_to_use, value = TRUE)
        if (length(ribo_genes) > 0) {
            genes_to_use <- setdiff(genes_to_use, ribo_genes)
            message("clusterCellTypes: excluded ", length(ribo_genes),
                    " ribosomal genes")
        }
    }

    # Subset to non-mito/ribo genes for normalization and clustering
    spe_sub <- spe[genes_to_use, ]

    # --- Normalize ---
    spe_sub <- scater::logNormCounts(spe_sub)

    # --- Covariate regression ---
    if (regress_covariates) {
        # Filter to covariates that actually exist in colData
        available_covs <- intersect(covariates, colnames(colData(spe)))
        missing_covs <- setdiff(covariates, colnames(colData(spe)))
        if (length(missing_covs) > 0) {
            message("clusterCellTypes: covariates not found, skipping: ",
                    paste(missing_covs, collapse = ", "))
        }

        if (length(available_covs) > 0) {
            message("clusterCellTypes: regressing out ",
                    paste(available_covs, collapse = ", "))

            # Build design matrix from covariates
            cov_df <- as.data.frame(colData(spe)[, available_covs,
                                                  drop = FALSE])
            # Scale covariates to avoid numerical issues
            cov_df <- as.data.frame(scale(cov_df))
            design <- model.matrix(~ ., data = cov_df)

            # Gene-by-gene regression on logcounts
            logcounts_mat <- as.matrix(assay(spe_sub, "logcounts"))

            # Compute residuals: Y - X %*% (X'X)^{-1} X' Y
            # Using lm.fit for each gene is simple and fast enough
            resid_mat <- apply(logcounts_mat, 1, function(y) {
                fit <- lm.fit(design, y)
                fit$residuals
            })
            resid_mat <- t(resid_mat)
            dimnames(resid_mat) <- dimnames(logcounts_mat)

            # Store regressed expression for PCA
            assay(spe_sub, "regressed") <- resid_mat
        }
    }

    # --- Feature selection ---
    use_assay <- if ("regressed" %in% SummarizedExperiment::assayNames(spe_sub)) {
        "regressed"
    } else {
        "logcounts"
    }

    dec <- scran::modelGeneVar(spe_sub, assay.type = use_assay)
    n_hvgs_actual <- min(n_hvgs, nrow(spe_sub))
    hvgs <- scran::getTopHVGs(dec, n = n_hvgs_actual)

    # --- PCA ---
    n_pcs_actual <- min(n_pcs, ncol(spe_sub) - 1, length(hvgs) - 1)
    spe_sub <- scater::runPCA(spe_sub, subset_row = hvgs,
                               ncomponents = n_pcs_actual,
                               exprs_values = use_assay)

    # --- Clustering ---
    clusters <- scran::clusterCells(spe_sub,
        use.dimred = "PCA",
        BLUSPARAM = bluster::NNGraphParam(
            cluster.fun = "leiden",
            cluster.args = list(resolution_parameter = resolution),
            k = k.neighbors
        )
    )

    cluster_labels <- as.character(clusters)

    # --- Handle small clusters ---
    cluster_sizes <- table(cluster_labels)
    small_clusters <- names(cluster_sizes[cluster_sizes < min_cluster_size])

    if (length(small_clusters) > 0) {
        n_reassigned <- sum(cluster_labels %in% small_clusters)
        cluster_labels[cluster_labels %in% small_clusters] <- "unassigned"
        message("clusterCellTypes: ", n_reassigned,
                " cells from ", length(small_clusters),
                " small cluster(s) labeled 'unassigned'")
    }

    # --- Store results ---
    colData(spe)[[cluster_col]] <- cluster_labels

    # Transfer PCA to the main object
    SingleCellExperiment::reducedDim(spe, "PCA") <-
        SingleCellExperiment::reducedDim(spe_sub, "PCA")

    n_clusters <- length(unique(cluster_labels))
    message("clusterCellTypes: identified ", n_clusters, " cluster(s)")

    spe
}
