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
#' For multi-sample datasets, batch effects are corrected using Harmony
#' (Korsunsky et al., Nature Methods 2019) on the PCA embedding before
#' clustering. This produces consistent cluster labels across samples
#' while downstream QC operates per-sample.
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
#' @param samples Column name in \code{colData} for sample IDs
#'   (default "sample_id"). When multiple samples are detected and the
#'   \code{harmony} package is installed, Harmony batch correction is applied
#'   to the PCA embedding before clustering. If the column is missing or
#'   only one sample is present, no batch correction is performed.
#' @param cluster_col Name of the \code{colData} column to store cluster
#'   labels (default "cell_cluster").
#' @param use_existing Name of an existing \code{colData} column containing
#'   pre-computed cluster labels. If provided, clustering is skipped entirely
#'   and these labels are used. Set to NULL (default) to perform clustering.
#' @param cluster_method Clustering algorithm to use. \code{"leiden"} (default)
#'   uses shared nearest-neighbor graph with Leiden community detection.
#'   \code{"mbkmeans"} uses mini-batch k-means, which is much faster for large
#'   datasets but requires specifying the number of clusters.
#' @param n_clusters Number of clusters for mini-batch k-means (default NULL).
#'   Only used when \code{cluster_method = "mbkmeans"}. If NULL, estimated
#'   automatically as \code{round(sqrt(ncol(spe) / 2))} capped between 5
#'   and 50.
#' @param subsample Maximum number of cells to cluster directly (default NULL).
#'   When set and \code{ncol(spe)} exceeds this value, a subsample is
#'   clustered and remaining cells are assigned to the nearest cluster centroid
#'   in PCA space. This "sketch" approach dramatically speeds up clustering for
#'   very large datasets (e.g., > 100k cells).
#' @param subsample_method Sampling method when \code{subsample} is used.
#'   \code{"leverage"} (default) samples cells with probability proportional
#'   to their PCA leverage scores, which preferentially selects cells in
#'   sparse regions of expression space (e.g., rare cell types).
#'   \code{"random"} uses uniform random sampling.
#' @param leverage_alpha Mixing weight between leverage and uniform sampling
#'   (default NULL = pure leverage). When set to a value in (0, 1), sampling
#'   probabilities are \code{alpha * leverage + (1 - alpha) * uniform}.
#'   Only used when \code{subsample_method = "leverage"}.
#' @param num_threads Number of threads for SNN graph construction and Leiden
#'   community detection (default 1). Passed to
#'   \code{bluster::NNGraphParam(num.threads = )}.
#' @param BPPARAM A \code{BiocParallelParam} object for parallel
#'   covariate regression and mini-batch k-means (default NULL = serial).
#' @param BSPARAM A \code{BiocSingularParam} object for PCA computation
#'   (default \code{BiocSingular::IrlbaParam()}, fast approximate PCA
#'   suitable for single-cell-resolution datasets). Set to
#'   \code{BiocSingular::ExactParam()} for an exact full SVD if results
#'   need to be bit-identical to base \code{prcomp}.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object for nearest-neighbor
#'   search in SNN graph construction (default NULL = exact kNN). Set to
#'   \code{BiocNeighbors::HnswParam()} for approximate NN.
#' @param seed Random seed for reproducibility (default 42).
#'
#' @return A \linkS4class{SpatialExperiment} object with cluster labels in
#'   \code{colData(spe)[[cluster_col]]}. PCA results are stored in
#'   \code{reducedDim(spe, "PCA")}. When Harmony batch correction is applied,
#'   the corrected embedding is stored in \code{reducedDim(spe, "HARMONY")}.
#'
#' @importFrom SummarizedExperiment colData colData<- assay assay<-
#' @importFrom SingleCellExperiment reducedDim reducedDim<-
#' @importFrom scran clusterCells modelGeneVar getTopHVGs
#' @importFrom scater logNormCounts runPCA
#' @importFrom bluster NNGraphParam MbkmeansParam
#' @importFrom BiocNeighbors queryKNN
#' @importFrom BiocSingular IrlbaParam ExactParam
#' @importFrom stats lm.fit model.matrix
#' @importFrom methods is
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
                             k.neighbors = 25,
                             samples = "sample_id",
                             cluster_col = "cell_cluster",
                             use_existing = NULL,
                             cluster_method = c("leiden", "mbkmeans"),
                             n_clusters = NULL,
                             subsample = NULL,
                             subsample_method = c("leverage", "random"),
                             leverage_alpha = NULL,
                             num_threads = 1,
                             BPPARAM = NULL,
                             BSPARAM = NULL,
                             BNPARAM = NULL,
                             seed = 42) {

    cluster_method <- match.arg(cluster_method)
    subsample_method <- match.arg(subsample_method)

    # --- Input validation ---
    # clusterCellTypes does not depend on spatial coordinates; accept
    # any SingleCellExperiment so users can reuse it on plain SCE data.
    if (!is(spe, "SingleCellExperiment")) {
        stop("'spe' must be a SingleCellExperiment or SpatialExperiment ",
            "object.")
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

    n_total <- ncol(spe)

    if (n_total < min_cluster_size) {
        warning("Fewer cells (", n_total, ") than min_cluster_size (",
            min_cluster_size, "). Assigning all cells to one cluster.")
        colData(spe)[[cluster_col]] <- "cluster_1"
        return(spe)
    }

    # Use withr to scope RNG state to this function call: the user's
    # RNG state is restored when the function returns. Inner re-seeds
    # at Harmony and sketch sites use the same idiom for the same
    # reason.
    withr::local_seed(seed)

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

    # Subset to filtered genes
    spe_work <- spe[genes_to_use, ]

    # --- Normalize ---
    spe_work <- scater::logNormCounts(spe_work)

    # --- Covariate regression (optional) ---
    if (regress_covariates) {
        available_covs <- intersect(covariates, colnames(colData(spe)))
        missing_covs <- setdiff(covariates, colnames(colData(spe)))
        if (length(missing_covs) > 0) {
            message("clusterCellTypes: covariates not found, skipping: ",
                paste(missing_covs, collapse = ", "))
        }
        if (length(available_covs) > 0) {
            message("clusterCellTypes: regressing out ",
                paste(available_covs, collapse = ", "))
            cov_df <- as.data.frame(
                colData(spe)[, available_covs, drop = FALSE])
            cov_df <- as.data.frame(scale(cov_df))
            design <- model.matrix(~., data = cov_df)

            # All-at-once regression: lm.fit accepts a multi-column
            # response matrix and solves the QR ONCE. This avoids the
            # gene-by-gene loop that would densify and replicate the
            # full logcounts matrix for every worker. Densification is
            # still required (lm.fit needs a numeric matrix), but it
            # happens exactly once and only the transposed view.
            lc <- assay(spe_work, "logcounts")
            Yt <- as.matrix(Matrix::t(lc))     # cells x genes
            fit <- stats::lm.fit(design, Yt)
            resid_mat <- t(fit$residuals)      # genes x cells
            dimnames(resid_mat) <- dimnames(lc)
            assay(spe_work, "regressed") <- resid_mat
            rm(Yt, fit)
            gc(verbose = FALSE)
        }
    }

    # --- Feature selection ---
    use_assay <- if ("regressed" %in%
        SummarizedExperiment::assayNames(spe_work)) {
        "regressed"
    } else {
        "logcounts"
    }

    dec <- scran::modelGeneVar(spe_work, assay.type = use_assay)
    n_hvgs_actual <- min(n_hvgs, nrow(spe_work))
    hvgs <- scran::getTopHVGs(dec, n = n_hvgs_actual)

    # --- PCA ---
    # `scater::runPCA` renamed `exprs_values` -> `assay.type`; the new
    # argument is the supported name on Bioc >= 3.18. Default to IRLBA
    # for tractable cost at single-cell-resolution scale.
    n_pcs_actual <- min(n_pcs, n_total - 1, length(hvgs) - 1)
    if (is.null(BSPARAM)) BSPARAM <- BiocSingular::IrlbaParam()
    spe_work <- scater::runPCA(spe_work,
        subset_row  = hvgs,
        ncomponents = n_pcs_actual,
        assay.type  = use_assay,
        BSPARAM     = BSPARAM)

    # Store PCA in original SPE
    SingleCellExperiment::reducedDim(spe, "PCA") <-
        SingleCellExperiment::reducedDim(spe_work, "PCA")

    # --- Harmony batch correction (multi-sample only) ---
    # Determine which reducedDim to use for clustering
    cluster_dimred <- "PCA"

    if (samples %in% colnames(colData(spe))) {
        sample_vec <- as.character(colData(spe)[[samples]])
        n_samples <- length(unique(sample_vec))
    } else {
        n_samples <- 1
    }

    if (n_samples > 1) {
        if (requireNamespace("harmony", quietly = TRUE)) {
            message("clusterCellTypes: running Harmony batch correction ",
                "across ", n_samples, " samples")
            pca_mat <- SingleCellExperiment::reducedDim(spe, "PCA")
            meta_df <- as.data.frame(colData(spe)[, samples, drop = FALSE])
            # Harmony has internal RNG with no seed kwarg; scope the seed
            # locally so user state is preserved.
            harmony_mat <- withr::with_seed(seed,
                harmony::RunHarmony(
                    data_mat  = pca_mat,
                    meta_data = meta_df,
                    vars_use  = samples,
                    verbose   = FALSE
                )
            )
            SingleCellExperiment::reducedDim(spe, "HARMONY") <- harmony_mat
            # Also store in spe_work so clusterCells can find it
            SingleCellExperiment::reducedDim(spe_work, "HARMONY") <-
                harmony_mat
            cluster_dimred <- "HARMONY"
        } else {
            message("clusterCellTypes: multiple samples detected but ",
                "'harmony' package not installed. Clustering without ",
                "batch correction. Install with: ",
                "install.packages('harmony')")
        }
    }

    # --- Subsampling (sketch) ---
    use_sub <- !is.null(subsample) && n_total > subsample
    if (use_sub) {
        pca_mat <- SingleCellExperiment::reducedDim(spe_work, cluster_dimred)

        # Make sketch index reproducible regardless of how much RNG
        # was consumed above (Harmony etc.), without polluting user RNG.
        sub_idx <- withr::with_seed(seed, {
            if (subsample_method == "leverage") {
                col_sd <- apply(pca_mat, 2, sd)
                col_sd[col_sd == 0] <- 1
                pca_scaled <- sweep(pca_mat, 2, col_sd, "/")
                leverage <- rowSums(pca_scaled^2)
                lev_prob <- leverage / sum(leverage)
                if (!is.null(leverage_alpha)) {
                    unif_prob <- rep(1 / n_total, n_total)
                    prob <- leverage_alpha * lev_prob +
                        (1 - leverage_alpha) * unif_prob
                } else {
                    prob <- lev_prob
                }
                message("clusterCellTypes: leverage-score subsampling ",
                    subsample, " of ", n_total, " cells")
                sample(n_total, subsample, prob = prob)
            } else {
                message("clusterCellTypes: random subsampling ",
                    subsample, " of ", n_total, " cells")
                sample(n_total, subsample)
            }
        })
        spe_sketch <- spe_work[, sub_idx]
    } else {
        spe_sketch <- spe_work
    }

    # --- Clustering ---
    if (cluster_method == "leiden") {
        nn_args <- list(cluster.fun = "leiden",
            cluster.args = list(
                resolution_parameter = resolution),
            k = min(k.neighbors, ncol(spe_sketch) - 1))
        if (num_threads > 1) nn_args$num.threads <- num_threads
        if (!is.null(BNPARAM)) nn_args$BNPARAM <- BNPARAM
        blusparam <- do.call(bluster::NNGraphParam, nn_args)
    } else {
        if (!requireNamespace("mbkmeans", quietly = TRUE)) {
            stop("Package 'mbkmeans' is required for ",
                "cluster_method = \"mbkmeans\". Install with:\n",
                "  BiocManager::install(\"mbkmeans\")")
        }
        n_cl_est <- n_clusters
        if (is.null(n_cl_est)) {
            n_cl_est <- max(5, min(50,
                round(sqrt(ncol(spe_sketch) / 2))))
            message("clusterCellTypes: n_clusters auto-estimated as ",
                n_cl_est)
        }
        mbk_args <- list(centers = n_cl_est)
        if (!is.null(BPPARAM) &&
            requireNamespace("BiocParallel", quietly = TRUE)) {
            mbk_args$BPPARAM <- BPPARAM
        }
        blusparam <- do.call(bluster::MbkmeansParam, mbk_args)
    }

    clusters <- scran::clusterCells(spe_sketch,
        use.dimred = cluster_dimred,
        BLUSPARAM = blusparam
    )
    cluster_labels <- as.character(clusters)

    # --- Project labels to remaining cells (sketch mode) ---
    if (use_sub) {
        emb_all <- SingleCellExperiment::reducedDim(spe_work, cluster_dimred)
        emb_sketch <- emb_all[sub_idx, , drop = FALSE]

        unique_cls <- unique(cluster_labels)
        centroids <- do.call(rbind, lapply(unique_cls, function(cl) {
            colMeans(emb_sketch[cluster_labels == cl, , drop = FALSE])
        }))

        nn_result <- BiocNeighbors::queryKNN(
            X = centroids, query = emb_all, k = 1
        )
        all_labels <- unique_cls[nn_result$index[, 1]]
    } else {
        all_labels <- cluster_labels
    }

    # --- Handle small clusters ---
    cluster_sizes <- table(all_labels)
    small_cls <- names(cluster_sizes[cluster_sizes < min_cluster_size])
    if (length(small_cls) > 0) {
        n_reassigned <- sum(all_labels %in% small_cls)
        all_labels[all_labels %in% small_cls] <- "unassigned"
        message("clusterCellTypes: ", n_reassigned, " cells from ",
            length(small_cls),
            " small cluster(s) labeled 'unassigned'")
    }

    # --- Store results ---
    colData(spe)[[cluster_col]] <- all_labels

    n_clusters_final <- length(unique(all_labels))
    message("clusterCellTypes: identified ", n_clusters_final,
        " cluster(s)")

    spe
}
