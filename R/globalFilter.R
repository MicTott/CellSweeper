#' Global Pre-Filtering of Low-Quality Cells
#'
#' Removes obvious debris and impossible segmentations from single-cell
#' resolution spatial transcriptomics data. This is Level 1 of the CellSweeper
#' QC framework and is intentionally permissive — it removes only clear junk
#' before clustering.
#'
#' @param spe A \linkS4class{SpatialExperiment} object.
#' @param min_counts Minimum total transcript counts per cell (default 5).
#'   Cells below this threshold are removed.
#' @param min_genes Minimum unique genes detected per cell (default 5).
#'   Cells below this threshold are removed.
#' @param max_area Maximum cell area. Cells above this threshold are removed.
#'   Set to NULL (default) to skip this filter.
#' @param min_area Minimum cell area. Cells below this threshold are removed.
#'   Set to NULL (default) to skip this filter.
#' @param max_neg_prop Maximum negative probe proportion. Cells above this
#'   threshold are removed. Set to NULL (default) to skip this filter.
#' @param counts_col Column name in \code{colData} for total counts
#'   (default "sum").
#' @param genes_col Column name in \code{colData} for unique genes detected
#'   (default "detected").
#' @param area_col Column name in \code{colData} for cell area
#'   (default "Area_um", matching SpaceTrooper output).
#' @param neg_col Column name in \code{colData} for negative probe proportion
#'   (default "altexps_NegPrb_percent", matching SpaceTrooper output).
#' @param add_transcript_density Logical. If TRUE (default) and \code{area_col}
#'   exists in colData, computes transcript density and adds \code{CountArea}
#'   and \code{log2CountArea} columns (matching SpaceTrooper naming). Skipped
#'   if \code{log2CountArea} already exists.
#'
#' @return A \linkS4class{SpatialExperiment} object with low-quality cells
#'   removed. A \code{globalFilter_pass} column is added to \code{colData}
#'   before subsetting.
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
#' counts <- matrix(rpois(1000, lambda = 5), nrow = 100, ncol = 10)
#' rownames(counts) <- paste0("Gene", 1:100)
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
#' spe <- globalFilter(spe, min_counts = 10, min_genes = 10)
globalFilter <- function(spe,
                         min_counts = 5,
                         min_genes = 5,
                         max_area = NULL,
                         min_area = NULL,
                         max_neg_prop = NULL,
                         counts_col = "sum",
                         genes_col = "detected",
                         area_col = "Area_um",
                         neg_col = "altexps_NegPrb_percent",
                         add_transcript_density = TRUE) {

    # --- Input validation ---
    if (!is(spe, "SpatialExperiment")) {
        stop("'spe' must be a SpatialExperiment object.")
    }

    cd <- colData(spe)
    cd_names <- colnames(cd)

    if (!counts_col %in% cd_names) {
        stop("Column '", counts_col, "' not found in colData. ",
             "Run scuttle::addPerCellQCMetrics() first.")
    }
    if (!genes_col %in% cd_names) {
        stop("Column '", genes_col, "' not found in colData. ",
             "Run scuttle::addPerCellQCMetrics() first.")
    }

    # --- Build keep mask ---
    n_start <- ncol(spe)
    keep <- rep(TRUE, n_start)

    # Minimum counts filter
    keep <- keep & (cd[[counts_col]] >= min_counts)

    # Minimum genes filter
    keep <- keep & (cd[[genes_col]] >= min_genes)

    # Area filters (only if specified AND column exists)
    if (!is.null(max_area) || !is.null(min_area)) {
        if (!area_col %in% cd_names) {
            message("Column '", area_col, "' not found in colData. ",
                    "Skipping area-based filtering.")
        } else {
            if (!is.null(max_area)) {
                keep <- keep & (cd[[area_col]] <= max_area)
            }
            if (!is.null(min_area)) {
                keep <- keep & (cd[[area_col]] >= min_area)
            }
        }
    }

    # Negative probe proportion filter
    if (!is.null(max_neg_prop) && !is.null(neg_col)) {
        if (!neg_col %in% cd_names) {
            message("Column '", neg_col, "' not found in colData. ",
                    "Skipping negative probe filtering.")
        } else {
            keep <- keep & (cd[[neg_col]] <= max_neg_prop)
        }
    }

    # --- Add transcript density (SpaceTrooper-compatible naming) ---
    if (add_transcript_density && area_col %in% cd_names) {
        if (!"log2CountArea" %in% cd_names) {
            raw_density <- cd[[counts_col]] / cd[[area_col]]
            colData(spe)$CountArea <- raw_density
            colData(spe)$log2CountArea <- log2(raw_density + 1)
        }
    }

    # --- Add filter flag and subset ---
    colData(spe)$globalFilter_pass <- keep

    n_removed <- sum(!keep)
    if (n_removed > 0) {
        message("globalFilter: removed ", n_removed, "/", n_start,
                " cells (", round(100 * n_removed / n_start, 1), "%)")
    } else {
        message("globalFilter: all ", n_start, " cells passed filters")
    }

    spe <- spe[, keep]
    spe
}
