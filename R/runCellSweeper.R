#' Run the Full CellSweeper QC Pipeline
#'
#' Convenience wrapper that runs all three levels of the CellSweeper QC
#' framework in sequence:
#'
#' \enumerate{
#'   \item \strong{Level 1:} Global pre-filtering (\code{\link{globalFilter}})
#'   \item \strong{Level 2a:} Clustering (\code{\link{clusterCellTypes}})
#'   \item \strong{Level 2b+c:} Artifact cluster detection
#'     (\code{\link{flagArtifactClusters}})
#'   \item \strong{Level 3:} Within-cluster spatial QC
#'     (\code{\link{clusterLocalOutliers}})
#' }
#'
#' @param spe A \linkS4class{SpatialExperiment} object with QC metrics in
#'   \code{colData} (e.g., from \code{scuttle::addPerCellQCMetrics()}).
#' @param min_counts,min_genes,min_area,max_area Thresholds passed to
#'   \code{\link{globalFilter}}.
#' @param resolution,regress_covariates,n_hvgs,use_existing,cluster_method,subsample,subsample_method,leverage_alpha,num_threads,BSPARAM
#'   Passed to \code{\link{clusterCellTypes}}.
#' @param n_permutations Number of permutations for the spatial coherence
#'   test (default 500). Passed to \code{\link{flagArtifactClusters}}.
#' @param metrics,morpho_metrics,k,z_threshold,combine,robust_z
#'   Passed to \code{\link{clusterLocalOutliers}}.
#' @param samples Column name for sample IDs (default \code{"sample_id"}).
#'   Shared across all stages.
#' @param BPPARAM,BNPARAM Shared parallelization and nearest-neighbor
#'   parameters; forwarded to clustering, artifact detection, and
#'   within-cluster outlier detection.
#' @param ... Additional arguments forwarded to
#'   \code{\link{clusterCellTypes}}, \code{\link{flagArtifactClusters}},
#'   and \code{\link{clusterLocalOutliers}}; useful for advanced options
#'   that the wrapper does not surface explicitly.
#' @param verbose Logical. Print progress messages (default TRUE).
#' @param seed Random seed for reproducibility (default 42).
#'
#' @return A \linkS4class{SpatialExperiment} object with all QC flags in
#'   \code{colData}. The final combined flag \code{cellsweeper_qc} is TRUE
#'   for cells that should be discarded (artifact cluster members OR
#'   within-cluster outliers).
#'
#' @importFrom SummarizedExperiment colData colData<-
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
#' spe$my_clusters <- rep("A", 10)
#'
#' # Run with pre-existing labels and minimal permutations
#' spe <- runCellSweeper(spe, use_existing = "my_clusters",
#'     morpho_metrics = NULL, n_permutations = 10, k = 5)
runCellSweeper <- function(spe,
                           # Level 1
                           min_counts = 5,
                           min_genes = 5,
                           max_area = NULL,
                           min_area = NULL,
                           # Level 2a
                           resolution = 0.5,
                           regress_covariates = TRUE,
                           n_hvgs = 2000,
                           use_existing = NULL,
                           cluster_method = c("leiden", "mbkmeans"),
                           subsample = NULL,
                           subsample_method = c("leverage", "random"),
                           leverage_alpha = NULL,
                           num_threads = 1,
                           # Level 2b+c
                           n_permutations = 500,
                           # Level 3
                           metrics = c("sum", "detected",
                               "subsets_mito_percent"),
                           morpho_metrics = c("Area_um",
                               "log2AspectRatio",
                               "log2CountArea"),
                           k = 50,
                           z_threshold = 3,
                           combine = c("union", "intersection"),
                           robust_z = c("mad", "iglewicz"),
                           # Common
                           samples = "sample_id",
                           BPPARAM = NULL,
                           BSPARAM = NULL,
                           BNPARAM = NULL,
                           verbose = TRUE,
                           seed = 42,
                           ...) {

    cluster_method   <- match.arg(cluster_method)
    subsample_method <- match.arg(subsample_method)
    combine          <- match.arg(combine)
    robust_z         <- match.arg(robust_z)

    if (!is(spe, "SpatialExperiment")) {
        stop("'spe' must be a SpatialExperiment object.")
    }

    n_start <- ncol(spe)

    # Partition `...` to the sub-call that accepts each named argument.
    dots <- list(...)
    .forward <- function(fn, extra = list()) {
        accepted <- names(formals(fn))
        keep <- intersect(names(dots), accepted)
        c(extra, dots[keep])
    }

    # --- Level 1: Global Pre-Filtering ---
    if (verbose) message("\n=== Level 1: Global pre-filtering ===")
    spe <- do.call(globalFilter, .forward(globalFilter, list(
        spe = spe,
        min_counts = min_counts, min_genes = min_genes,
        max_area = max_area, min_area = min_area)))

    # --- Level 2a: Clustering ---
    if (verbose) message("\n=== Level 2a: Clustering cell types ===")
    spe <- do.call(clusterCellTypes, .forward(clusterCellTypes, list(
        spe = spe,
        resolution = resolution,
        regress_covariates = regress_covariates,
        n_hvgs = n_hvgs,
        samples = samples,
        use_existing = use_existing,
        cluster_method = cluster_method,
        subsample = subsample,
        subsample_method = subsample_method,
        leverage_alpha = leverage_alpha,
        num_threads = num_threads,
        BPPARAM = BPPARAM, BSPARAM = BSPARAM, BNPARAM = BNPARAM,
        seed = seed)))

    # --- Level 2b+c: Flag Artifact Clusters ---
    if (verbose) message("\n=== Level 2b: Flagging artifact clusters ===")
    spe <- do.call(flagArtifactClusters, .forward(flagArtifactClusters, list(
        spe = spe,
        metrics = metrics,
        n_permutations = n_permutations,
        samples = samples,
        BPPARAM = BPPARAM, BNPARAM = BNPARAM,
        seed = seed)))

    # --- Level 3: Within-Cluster Spatial QC ---
    if (verbose) message("\n=== Level 3: Within-cluster spatial QC ===")
    spe <- do.call(clusterLocalOutliers, .forward(clusterLocalOutliers, list(
        spe = spe,
        metrics = metrics,
        morpho_metrics = morpho_metrics,
        k = k,
        z_threshold = z_threshold,
        combine = combine,
        robust_z = robust_z,
        samples = samples,
        BPPARAM = BPPARAM, BNPARAM = BNPARAM)))

    # --- Final combined flag ---
    cd <- colData(spe)
    cd$cellsweeper_qc <- cd$cluster_artifact_flag | cd$cellsweeper_outlier
    # Handle NAs
    cd$cellsweeper_qc[is.na(cd$cellsweeper_qc)] <- FALSE
    colData(spe) <- cd

    if (verbose) {
        n_final <- ncol(spe)
        n_flagged <- sum(cd$cellsweeper_qc, na.rm = TRUE)
        n_removed_global <- n_start - n_final
        message("\n=== CellSweeper Summary ===")
        message("  Cells input:           ", n_start)
        message("  Removed (global):      ", n_removed_global)
        message("  Flagged (cluster QC):  ",
            sum(cd$cluster_artifact_flag, na.rm = TRUE))
        message("  Flagged (spatial QC):  ",
            sum(cd$cellsweeper_outlier, na.rm = TRUE))
        message("  Total flagged:         ", n_flagged, "/", n_final,
            " (", round(100 * n_flagged / n_final, 1), "%)")
    }

    spe
}
