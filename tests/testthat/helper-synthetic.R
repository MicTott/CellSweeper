# Shared test data factory for CellSweeper tests
# Creates a synthetic SpatialExperiment with known structure

#' Create a synthetic SpatialExperiment for testing
#'
#' Builds a small SPE with 3 spatial clusters (blobs in 2D),
#' known QC metrics, and optionally injected outliers.
#'
#' @param n_cells_per_cluster Number of cells per cluster (default 150)
#' @param n_genes Number of genes (default 200)
#' @param n_outliers Number of outlier cells to inject per cluster (default 5)
#' @param include_morphology Whether to add morphology columns (default TRUE)
#' @param include_polygons Whether to add synthetic sf polygons (default FALSE)
#' @param seed Random seed for reproducibility
#' @return A SpatialExperiment object
make_test_spe <- function(n_cells_per_cluster = 150,
                          n_genes = 200,
                          n_outliers = 5,
                          include_morphology = TRUE,
                          include_polygons = FALSE,
                          seed = 42) {

    set.seed(seed)

    n_clusters <- 3
    n_cells <- n_cells_per_cluster * n_clusters

    # --- Spatial coordinates: 3 blobs ---
    centers <- matrix(c(0, 0, 100, 0, 50, 100), ncol = 2, byrow = TRUE)
    coords <- do.call(rbind, lapply(seq_len(n_clusters), function(i) {
        matrix(rnorm(n_cells_per_cluster * 2, mean = 0, sd = 10),
            ncol = 2) +
            matrix(rep(centers[i, ], each = n_cells_per_cluster), ncol = 2)
    }))
    colnames(coords) <- c("x_centroid", "y_centroid")

    # --- Cluster labels ---
    cluster_labels <- rep(paste0("cluster_", seq_len(n_clusters)),
        each = n_cells_per_cluster)

    # --- Gene expression ---
    base_means <- c(5, 3, 1)
    counts_matrix <- do.call(cbind, lapply(seq_len(n_clusters), function(i) {
        matrix(rpois(n_genes * n_cells_per_cluster,
            lambda = base_means[i]),
        nrow = n_genes, ncol = n_cells_per_cluster)
    }))

    # Cluster-specific marker genes (first 10 per cluster)
    for (i in seq_len(n_clusters)) {
        gene_start <- (i - 1) * 10 + 1
        gene_end <- i * 10
        col_start <- (i - 1) * n_cells_per_cluster + 1
        col_end <- i * n_cells_per_cluster
        counts_matrix[gene_start:gene_end, col_start:col_end] <-
            counts_matrix[gene_start:gene_end, col_start:col_end] + 10
    }

    # Mito genes: last 10
    mito_idx <- (n_genes - 9):n_genes
    gene_names <- paste0("Gene", seq_len(n_genes))
    gene_names[mito_idx] <- paste0("MT-", gene_names[mito_idx])
    rownames(counts_matrix) <- gene_names
    colnames(counts_matrix) <- paste0("cell_", seq_len(n_cells))

    # --- Inject outliers INTO the counts matrix BEFORE creating SPE ---
    outlier_flags <- rep(FALSE, n_cells)
    if (n_outliers > 0) {
        for (i in seq_len(n_clusters)) {
            cluster_start <- (i - 1) * n_cells_per_cluster + 1
            outlier_idx <- cluster_start:(cluster_start + n_outliers - 1)
            outlier_flags[outlier_idx] <- TRUE

            # Near-zero library size
            counts_matrix[, outlier_idx] <- 0
            counts_matrix[1:2, outlier_idx] <- 1

            # High mito in a different cell
            mid_outlier <- cluster_start + n_outliers
            if (mid_outlier <= i * n_cells_per_cluster) {
                counts_matrix[mito_idx, mid_outlier] <-
                    counts_matrix[mito_idx, mid_outlier] * 10
            }
        }
    }

    # --- Build SpatialExperiment (with modified counts) ---
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = as(counts_matrix, "dgCMatrix")),
        spatialCoords = coords,
        sample_id = rep("sample1", n_cells)
    )

    # --- Compute QC metrics ONCE ---
    spe <- scuttle::addPerCellQCMetrics(spe,
        subsets = list(mito = grep("^MT-", rownames(spe))))

    # --- Add morphology columns (matching SpaceTrooper output names) ---
    if (include_morphology) {
        SummarizedExperiment::colData(spe)$Area_um <-
            pmax(50 + spe$sum * 0.05 + rnorm(n_cells, 0, 10), 10)
        SummarizedExperiment::colData(spe)$log2AspectRatio <-
            rnorm(n_cells, 0, 0.5)
        raw_density <- spe$sum /
            SummarizedExperiment::colData(spe)$Area_um
        SummarizedExperiment::colData(spe)$CountArea <- raw_density
        SummarizedExperiment::colData(spe)$log2CountArea <-
            log2(raw_density + 1)

        # Morphology outliers
        if (n_outliers > 0) {
            for (i in seq_len(n_clusters)) {
                cluster_start <- (i - 1) * n_cells_per_cluster + 1
                outlier_idx <- cluster_start:(cluster_start + n_outliers - 1)
                SummarizedExperiment::colData(spe)$Area_um[outlier_idx] <- 5
                SummarizedExperiment::colData(spe)$log2AspectRatio[outlier_idx] <- 3.0
            }
        }
    }

    # --- Store ground truth ---
    SummarizedExperiment::colData(spe)$true_outlier <- outlier_flags
    SummarizedExperiment::colData(spe)$true_cluster <- cluster_labels

    # --- Optional: synthetic polygons for plot testing ---
    if (include_polygons && requireNamespace("sf", quietly = TRUE)) {
        sp_coords <- SpatialExperiment::spatialCoords(spe)
        poly_list <- lapply(seq_len(nrow(sp_coords)), function(i) {
            x <- sp_coords[i, 1]
            y <- sp_coords[i, 2]
            # Small square around each centroid
            sf::st_polygon(list(matrix(c(
                x - 1, y - 1,
                x + 1, y - 1,
                x + 1, y + 1,
                x - 1, y + 1,
                x - 1, y - 1
            ), ncol = 2, byrow = TRUE)))
        })
        SummarizedExperiment::colData(spe)$polygons <-
            sf::st_sfc(poly_list)
    }

    spe
}


#' Create a multi-sample synthetic SpatialExperiment for testing
#'
#' Builds a 2-sample SPE by generating data for two samples with different
#' seeds and spatial offsets, then constructing a single SPE.
#'
#' @param n_cells_per_cluster Number of cells per cluster per sample (default 100)
#' @param n_genes Number of genes (default 200)
#' @param n_outliers Number of outlier cells per cluster per sample (default 3)
#' @return A SpatialExperiment with sample_id = "sample1" / "sample2"
make_test_spe_multisample <- function(n_cells_per_cluster = 100,
                                      n_genes = 200,
                                      n_outliers = 3) {
    # Build each sample's raw data independently
    .build_sample_data <- function(seed, spatial_x_offset = 0) {
        set.seed(seed)
        n_clusters <- 3
        n_cells <- n_cells_per_cluster * n_clusters

        # Spatial coordinates: 3 blobs
        centers <- matrix(c(0, 0, 100, 0, 50, 100), ncol = 2, byrow = TRUE)
        coords <- do.call(rbind, lapply(seq_len(n_clusters), function(i) {
            matrix(rnorm(n_cells_per_cluster * 2, mean = 0, sd = 10),
                ncol = 2) +
                matrix(rep(centers[i, ], each = n_cells_per_cluster), ncol = 2)
        }))
        coords[, 1] <- coords[, 1] + spatial_x_offset
        colnames(coords) <- c("x_centroid", "y_centroid")

        # Cluster labels
        cluster_labels <- rep(paste0("cluster_", seq_len(n_clusters)),
            each = n_cells_per_cluster)

        # Gene expression
        base_means <- c(5, 3, 1)
        counts_matrix <- do.call(cbind, lapply(seq_len(n_clusters), function(i) {
            matrix(rpois(n_genes * n_cells_per_cluster,
                lambda = base_means[i]),
            nrow = n_genes, ncol = n_cells_per_cluster)
        }))
        for (i in seq_len(n_clusters)) {
            gene_start <- (i - 1) * 10 + 1
            gene_end <- i * 10
            col_start <- (i - 1) * n_cells_per_cluster + 1
            col_end <- i * n_cells_per_cluster
            counts_matrix[gene_start:gene_end, col_start:col_end] <-
                counts_matrix[gene_start:gene_end, col_start:col_end] + 10
        }

        # Mito genes
        mito_idx <- (n_genes - 9):n_genes
        gene_names <- paste0("Gene", seq_len(n_genes))
        gene_names[mito_idx] <- paste0("MT-", gene_names[mito_idx])
        rownames(counts_matrix) <- gene_names

        # Inject outliers
        outlier_flags <- rep(FALSE, n_cells)
        if (n_outliers > 0) {
            for (i in seq_len(n_clusters)) {
                cluster_start <- (i - 1) * n_cells_per_cluster + 1
                outlier_idx <- cluster_start:(cluster_start + n_outliers - 1)
                outlier_flags[outlier_idx] <- TRUE
                counts_matrix[, outlier_idx] <- 0
                counts_matrix[1:2, outlier_idx] <- 1
                mid_outlier <- cluster_start + n_outliers
                if (mid_outlier <= i * n_cells_per_cluster) {
                    counts_matrix[mito_idx, mid_outlier] <-
                        counts_matrix[mito_idx, mid_outlier] * 10
                }
            }
        }

        list(counts = counts_matrix, coords = coords,
            cluster_labels = cluster_labels, outlier_flags = outlier_flags,
            n_cells = n_cells)
    }

    d1 <- .build_sample_data(seed = 42, spatial_x_offset = 0)
    d2 <- .build_sample_data(seed = 123, spatial_x_offset = 500)

    n_total <- d1$n_cells + d2$n_cells

    # Combine counts
    counts_all <- cbind(d1$counts, d2$counts)
    colnames(counts_all) <- c(paste0("cell_", seq_len(d1$n_cells)),
        paste0("s2_cell_", seq_len(d2$n_cells)))

    # Combine coordinates
    coords_all <- rbind(d1$coords, d2$coords)

    # Sample IDs
    sample_ids <- c(rep("sample1", d1$n_cells), rep("sample2", d2$n_cells))

    # Build combined SPE
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = as(counts_all, "dgCMatrix")),
        spatialCoords = coords_all,
        sample_id = sample_ids
    )

    # Compute QC metrics
    spe <- scuttle::addPerCellQCMetrics(spe,
        subsets = list(mito = grep("^MT-", rownames(spe))))

    # Add morphology columns
    n_cells_vec <- c(d1$n_cells, d2$n_cells)
    SummarizedExperiment::colData(spe)$Area_um <-
        pmax(50 + spe$sum * 0.05 + rnorm(n_total, 0, 10), 10)
    SummarizedExperiment::colData(spe)$log2AspectRatio <-
        rnorm(n_total, 0, 0.5)
    raw_density <- spe$sum / SummarizedExperiment::colData(spe)$Area_um
    SummarizedExperiment::colData(spe)$CountArea <- raw_density
    SummarizedExperiment::colData(spe)$log2CountArea <- log2(raw_density + 1)

    # Morphology outliers
    if (n_outliers > 0) {
        for (sample_offset in c(0, d1$n_cells)) {
            for (i in seq_len(3)) {
                cluster_start <- sample_offset + (i - 1) * n_cells_per_cluster + 1
                outlier_idx <- cluster_start:(cluster_start + n_outliers - 1)
                SummarizedExperiment::colData(spe)$Area_um[outlier_idx] <- 5
                SummarizedExperiment::colData(spe)$log2AspectRatio[outlier_idx] <- 3.0
            }
        }
    }

    # Ground truth
    SummarizedExperiment::colData(spe)$true_outlier <-
        c(d1$outlier_flags, d2$outlier_flags)
    SummarizedExperiment::colData(spe)$true_cluster <-
        c(d1$cluster_labels, d2$cluster_labels)

    spe
}
