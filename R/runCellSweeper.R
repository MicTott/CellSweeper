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
#' @param min_counts Minimum total counts for global filter (default 5).
#' @param min_genes Minimum unique genes for global filter (default 5).
#' @param max_area Maximum cell area for global filter (default NULL).
#' @param min_area Minimum cell area for global filter (default NULL).
#' @param resolution Leiden clustering resolution (default 0.5).
#' @param regress_covariates Logical. Regress QC covariates before PCA
#'   (default TRUE).
#' @param n_hvgs Number of HVGs for clustering (default 2000).
#' @param use_existing Name of existing \code{colData} column with cluster
#'   labels (default NULL, perform clustering).
#' @param metrics Transcriptomic QC metrics for outlier detection
#'   (default: \code{c("sum", "detected", "subsets_mito_percent")}).
#' @param morpho_metrics Morphological metrics for outlier detection
#'   (default: \code{c("Area_um", "log2AspectRatio", "log2CountArea")},
#'   matching SpaceTrooper output). Missing columns are skipped.
#' @param k Number of spatial nearest neighbors for within-cluster QC
#'   (default 50).
#' @param z_threshold Z-score threshold for outlier detection (default 3).
#' @param n_permutations Number of permutations for spatial coherence test
#'   (default 500).
#' @param combine Method for combining outlier flags: \code{"union"}
#'   (default) or \code{"intersection"}.
#' @param samples Column name for sample IDs (default "sample_id").
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
#'                       morpho_metrics = NULL, n_permutations = 10, k = 5)
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
                           combine = "union",
                           # Common
                           samples = "sample_id",
                           verbose = TRUE,
                           seed = 42) {

    if (!is(spe, "SpatialExperiment")) {
        stop("'spe' must be a SpatialExperiment object.")
    }

    n_start <- ncol(spe)

    # --- Level 1: Global Pre-Filtering ---
    if (verbose) message("\n=== Level 1: Global pre-filtering ===")
    spe <- globalFilter(spe,
                        min_counts = min_counts,
                        min_genes = min_genes,
                        max_area = max_area,
                        min_area = min_area)

    # --- Level 2a: Clustering ---
    if (verbose) message("\n=== Level 2a: Clustering cell types ===")
    spe <- clusterCellTypes(spe,
                            resolution = resolution,
                            regress_covariates = regress_covariates,
                            n_hvgs = n_hvgs,
                            use_existing = use_existing,
                            seed = seed)

    # --- Level 2b+c: Flag Artifact Clusters ---
    if (verbose) message("\n=== Level 2b: Flagging artifact clusters ===")
    spe <- flagArtifactClusters(spe,
                                metrics = metrics,
                                n_permutations = n_permutations,
                                samples = samples,
                                seed = seed)

    # --- Level 3: Within-Cluster Spatial QC ---
    if (verbose) message("\n=== Level 3: Within-cluster spatial QC ===")
    spe <- clusterLocalOutliers(spe,
                                metrics = metrics,
                                morpho_metrics = morpho_metrics,
                                k = k,
                                z_threshold = z_threshold,
                                combine = combine,
                                samples = samples)

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
