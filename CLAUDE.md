# CellSweeper: Cluster-Aware Spatial Quality Control for Single-Cell Resolution Spatial Transcriptomics

## Project Overview

CellSweeper is a new R/Bioconductor package for quality control of **single-cell resolution spatial transcriptomics** data (Xenium, CosMx, MERFISH, Stereo-seq, Slide-seq, VisiumHD). It extends the SpotSweeper framework (Totty, Hicks & Guo, Nature Methods 2025) to single-cell resolution, where the core assumption of SpotSweeper — that spatial neighbors are biologically similar — breaks down because adjacent cells can be entirely different types (e.g., neurons next to glia).

**Core innovation:** CellSweeper clusters cells in gene expression space to group biologically similar cells, then performs SpotSweeper-style spatial outlier detection *within each cluster* in physical tissue coordinates. This restores the local homogeneity assumption required for robust outlier detection. The method additionally incorporates segmentation morphology metrics and a cluster-level QC step to identify entire artifact clusters.

## Motivation & Background

### Why SpotSweeper doesn't directly apply to single-cell resolution data
- SpotSweeper computes local z-scores on QC metrics (library size, unique genes, mito %) using k-nearest spatial neighbors
- This works for Visium because ~55µm spots sample relatively homogeneous local biology — spatial neighbors tend to be in the same tissue domain
- At single-cell resolution, a neuron (high library size) can sit directly adjacent to a glial cell (lower library size) — spatial neighbors are biologically heterogeneous
- Comparing QC metrics across cell types inflates false positives (biological variation mistaken for technical outliers) and misses true outliers within cell types

### Why existing cluster-specific QC methods are insufficient
- **ddqc** (Subramanian et al., Genome Biology 2022): Clusters in expression space, applies per-cluster MAD thresholds. Works for scRNA-seq but has NO spatial component — cannot detect spatially-correlated artifacts or leverage spatial context
- **ctQC** (bioRxiv 2024): Similar cell-type-specific approach, showed spatial coherence improves post-QC on Slide-seq, but does not deeply integrate spatial or morphological features
- **Neither method addresses:** segmentation quality, regional spatial artifacts in imaging-based platforms, or the problem of entire low-quality clusters

### CellSweeper's approach
CellSweeper fills this gap with a **three-level QC framework**:

1. **Global pre-filtering** — Remove obvious debris and impossible segmentations
2. **Cluster-level QC** — Identify entire artifact clusters via multivariate pseudobulk QC outlier detection + spatial dispersion analysis
3. **Within-cluster spatial QC** — SpotSweeper-style local outlier detection in tissue space, restricted to same-cluster cells, on both transcriptomic and morphological metrics

## Method Details

### Level 1: Global Pre-Filtering
- Remove cells with near-zero transcript counts (e.g., < 5-10 transcripts)
- Remove cells with physically impossible segmentations (extremely small/large area)
- Platform-specific: high negative probe proportion (Xenium, CosMx)
- This is intentionally permissive — just removing obvious junk before clustering

### Level 2: Cluster-Level QC (Identifying Artifact Clusters)

#### 2a. Clustering
- Perform coarse, low-resolution clustering in gene expression space
- Goal: broad cell-type groupings (neurons vs glia vs endothelial), NOT fine subtypes
- Use standard Bioconductor workflow: normalize → HVGs → PCA → shared nearest-neighbor graph → Leiden/Louvain clustering (via `bluster`)
- **Important:** QC metrics (library size, mito %, etc.) can contaminate PCA. Consider:
  - Excluding mito/ribo genes from feature selection before PCA
  - Residualizing gene expression on QC covariates (library size, cell area) before PCA — regress out QC metrics at the gene level, THEN do PCA and clustering. This is analogous to `ScaleData(vars.to.regress = ...)` in Seurat
  - Using low clustering resolution to get robust, coarse groupings
- Minimum cluster size threshold (e.g., 50-100 cells) — clusters smaller than this fall back to global QC
- Alternatively, support user-provided cluster labels or reference-based label transfer (e.g., from SingleR) as input

#### 2b. Pseudobulk Multivariate Outlier Detection
- For each cluster, compute summary statistics: median library size, median unique genes, median mito %, median cell area, median eccentricity, median transcript density (counts/area), median negative probe proportion
- Represent each cluster as a feature vector of these summary statistics
- Compute **robust Mahalanobis distance** of each cluster from the multivariate center (use minimum covariance determinant estimator for robustness with few clusters)
- Flag clusters with high Mahalanobis distance as potential artifact clusters
- Also consider within-cluster variance of QC metrics as additional signal

#### 2c. Spatial Dispersion Analysis
- Real cell types are spatially organized in tissue; artifact clusters (driven by shared damage signatures) tend to be spatially scattered
- For each cluster, compute **neighborhood homogeneity**: for each cell, what proportion of its k spatial nearest neighbors share its cluster label? Average across the cluster
- Compare observed neighborhood homogeneity against a **permutation null** (randomly permute cluster labels) to account for cluster size differences — large clusters will have high homogeneity by chance
- Compute a z-score or p-value for each cluster's spatial coherence relative to null
- A cluster that is BOTH a multivariate QC outlier AND spatially incoherent is almost certainly an artifact
- A cluster with low library size but strong spatial coherence may be a legitimate low-expression cell type (e.g., quiescent cells in a niche) — do NOT flag these

### Level 3: Within-Cluster Spatial QC (Cell-Level Outlier Detection)

#### 3a. Transcriptomic QC Metrics
- For each cluster that passed Level 2, perform SpotSweeper-style local outlier detection:
  - Define k-nearest neighbors in **physical tissue coordinates** (Euclidean distance), restricted to cells within the same cluster
  - Compute local z-scores for: library size, unique genes detected, mitochondrial percentage
  - Flag cells with z-scores beyond threshold (default: |z| > 3) as outliers
- This is the direct generalization of SpotSweeper: spatial neighborhoods are now biologically homogeneous because we pre-filtered by cluster

#### 3b. Morphological/Segmentation QC Metrics
- Apply the same within-cluster spatial outlier detection on segmentation metrics:
  - **Cell/nucleus area** — too small = fragmented, too large = merged cells
  - **Eccentricity** — unusually elongated shapes suggest segmentation errors
  - **Compactness** (area / perimeter^2) — irregular shapes are suspect
  - **Transcript density** (counts per unit area) — partially decouples library size from segmentation quality
  - **Negative probe proportion** (platform-specific) — non-specific binding indicator
  - **Distance to FOV borders** (CosMx especially) — stitching artifacts at field-of-view boundaries
- Some morphological metrics may be less cell-type-dependent (a badly segmented cell looks weird regardless of type), so consider also offering a global morphological QC pass in addition to within-cluster

#### 3c. Combined Outlier Flag
- A cell is flagged as low quality if it is an outlier in ANY metric (union approach, as in SpotSweeper) or optionally require outlier status in multiple metrics (intersection approach)
- Store per-metric outlier flags and a combined flag in colData for user inspection

## Package Architecture

### Bioconductor Framework
- **Extends SpatialExperiment** (SPE) — inherits from SingleCellExperiment, includes spatial coordinates
- All results stored in `colData()` of the SPE object, consistent with SpotSweeper
- Uses existing Bioconductor infrastructure:
  - `BiocNeighbors` for k-nearest neighbor detection (both spatial and expression space)
  - `scuttle` for `addPerCellQCMetrics()` — standard QC metric computation
  - `scran` for normalization (`computePooledFactors`, `logNormCounts`)
  - `scater` for PCA (`runPCA`)
  - `bluster` for clustering (`clusterCells`, `NNGraphParam`)
  - `robustbase` or similar for robust Mahalanobis distance
  - `escheR` for spatial visualization (as in SpotSweeper)

### Naming Conventions (following SpotSweeper patterns)
- SpotSweeper functions: `localOutliers()`, `localVariance()`, `findArtifacts()`, `flagVisiumOutliers()`, `plotQCmetrics()`, `plotQCpdf()`
- CellSweeper should follow similar patterns with clear, descriptive function names

### Key Functions to Implement

```r
# ============================================================
# Level 1: Global Pre-Filtering
# ============================================================

globalFilter(spe,
    min_counts = 5,          # minimum total transcript counts
    min_genes = 5,           # minimum unique genes detected
    max_area = NULL,         # maximum cell area (platform-dependent)
    min_area = NULL,         # minimum cell area
    max_neg_prop = NULL,     # maximum negative probe proportion
    area_col = "cell_area",  # colData column for cell area
    neg_col = NULL           # colData column for negative probe counts
)
# Returns: SPE with obvious junk removed + colData flag

# ============================================================
# Level 2: Cluster-Level QC
# ============================================================

clusterCellTypes(spe,
    resolution = 0.5,            # Leiden resolution (low for coarse types)
    regress_covariates = TRUE,   # regress QC metrics from expression before PCA
    covariates = c("sum", "detected", "subsets_mito_percent", "cell_area"),
    n_hvgs = 2000,
    n_pcs = 30,
    min_cluster_size = 50,       # minimum cells per cluster
    cluster_col = "cell_cluster", # output colData column
    use_existing = NULL          # optionally use pre-existing cluster labels
)
# Returns: SPE with cluster labels in colData
# If use_existing is provided, skip clustering and use those labels

flagArtifactClusters(spe,
    cluster_col = "cell_cluster",
    metrics = c("sum", "detected", "subsets_mito_percent",
                "cell_area", "transcript_density"),
    # Pseudobulk multivariate outlier detection
    mahal_threshold = "auto",     # or numeric chi-squared threshold
    covariance_method = "mcd",    # minimum covariance determinant
    # Spatial dispersion
    spatial_k = 20,               # k neighbors for spatial homogeneity
    n_permutations = 500,         # permutations for null distribution
    dispersion_threshold = 0.05,  # p-value threshold for spatial incoherence
    # Combined: flag if BOTH QC outlier AND spatially incoherent
    require_both = TRUE
)
# Returns: SPE with cluster-level flags in colData
# Adds: "cluster_qc_mahal", "cluster_spatial_homogeneity",
#        "cluster_spatial_pvalue", "cluster_artifact_flag"

# ============================================================
# Level 3: Within-Cluster Spatial QC
# ============================================================

clusterLocalOutliers(spe,
    cluster_col = "cell_cluster",
    exclude_artifact_clusters = TRUE,  # skip clusters flagged in Level 2
    # Transcriptomic metrics
    metrics = c("sum", "detected", "subsets_mito_percent"),
    # Morphological metrics (NULL = skip if not available)
    morpho_metrics = c("cell_area", "eccentricity", "compactness",
                       "transcript_density"),
    # Spatial neighborhood parameters
    k = 50,                    # spatial nearest neighbors (within cluster)
    min_neighbors = 10,        # minimum same-cluster neighbors required
    # Outlier detection
    z_threshold = 3,           # z-score threshold for outlier detection
    direction = "both",        # "lower", "upper", or "both" per metric
    log_transform = TRUE,      # log-transform skewed metrics before z-score
    # Output
    combine = "union"          # "union" or "intersection" for multi-metric
)
# Returns: SPE with per-metric outlier flags and combined flag in colData
# Adds: "{metric}_zscore", "{metric}_outlier", "cellsweeper_outlier"

# ============================================================
# Convenience / Wrapper
# ============================================================

runCellSweeper(spe,
    # Runs the full pipeline: globalFilter → clusterCellTypes →
    # flagArtifactClusters → clusterLocalOutliers
    # Accepts all parameters from individual functions
    ...
)
# Returns: fully QC'd SPE with all flags in colData

# ============================================================
# Visualization
# ============================================================

plotClusterQC(spe,
    cluster_col = "cell_cluster",
    metric = "sum",
    # Shows per-cluster distributions with outlier thresholds
)

plotSpatialOutliers(spe,
    metric = "sum",
    outlier_col = "cellsweeper_outlier",
    # Spatial scatter plot highlighting outlier cells
)

plotSpatialDispersion(spe,
    cluster_col = "cell_cluster",
    # Spatial plot colored by cluster, highlighting artifact clusters
)

plotClusterSummary(spe,
    cluster_col = "cell_cluster",
    # Panel showing: cluster-level QC metrics, spatial homogeneity scores,
    # Mahalanobis distances, artifact flags
)
```

### colData Output Convention
Following SpotSweeper, all results go into `colData(spe)`:
- `cell_cluster` — cluster assignment from Level 2
- `cluster_artifact_flag` — TRUE/FALSE per cell (inherited from cluster-level QC)
- `cluster_qc_mahal` — Mahalanobis distance of cell's cluster
- `cluster_spatial_homogeneity` — neighborhood homogeneity score of cell's cluster
- `cluster_spatial_pvalue` — permutation p-value for spatial coherence
- `{metric}_zscore` — local z-score for each QC metric (Level 3)
- `{metric}_outlier` — TRUE/FALSE outlier flag per metric
- `cellsweeper_outlier` — combined outlier flag (union/intersection of all metrics)
- `cellsweeper_qc` — final combined flag (TRUE = discard) incorporating all three levels

## Package Structure

```
CellSweeper/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── globalFilter.R           # Level 1: global pre-filtering
│   ├── clusterCellTypes.R       # Level 2a: expression clustering with QC regression
│   ├── flagArtifactClusters.R   # Level 2b+2c: pseudobulk outlier + spatial dispersion
│   ├── clusterLocalOutliers.R   # Level 3: within-cluster spatial outlier detection
│   ├── runCellSweeper.R         # Convenience wrapper for full pipeline
│   ├── spatialHomogeneity.R     # Spatial dispersion/homogeneity computation
│   ├── utils.R                  # Shared utility functions
│   ├── plot-functions.R         # Visualization functions
│   └── AllGenerics.R            # S4 method definitions if needed
├── man/                         # roxygen2-generated documentation
├── vignettes/
│   └── getting_started.Rmd      # BiocStyle vignette
├── tests/
│   └── testthat/                # Unit tests
├── inst/
│   └── extdata/                 # Example datasets if needed
└── .github/
    └── workflows/               # CI via GitHub Actions (BiocCheck, R CMD check)
```

## Dependencies

### Imports (required)
- `SpatialExperiment` — core data structure
- `SingleCellExperiment` — parent class
- `SummarizedExperiment` — base class
- `S4Vectors` — DataFrame, colData
- `BiocNeighbors` — k-NN computation (both spatial and expression)
- `scuttle` — `addPerCellQCMetrics`, `isOutlier`, `perCellQCFilters`
- `scran` — `computePooledFactors`, `modelGeneVar`, `getTopHVGs`
- `scater` — `runPCA`, `logNormCounts`
- `bluster` — `clusterCells`, `NNGraphParam`
- `stats` — basic statistical functions
- `Matrix` — sparse matrix operations

### Suggests (optional)
- `SpaceTrooper` — recommended for computing morphology/segmentation QC metrics from platform-specific data (Xenium, CosMx, VisiumHD, etc.). Users should preprocess with SpaceTrooper before running CellSweeper
- `robustbase` — `covMcd` for robust covariance estimation (Mahalanobis distance)
- `escheR` — spatial visualization
- `ggplot2` — plotting
- `BiocStyle` — vignette styling
- `testthat` — unit testing
- `STexampleData` — example datasets for vignettes
- `knitr`, `rmarkdown` — vignette building

## Implementation Notes

### Morphology Metrics and SpaceTrooper

CellSweeper **does NOT handle I/O or computation of morphology/segmentation metrics** (cell area, eccentricity, compactness, FOV border distances, etc.). These are non-trivial to extract and are highly platform-specific (Xenium HDF5, CosMx FOV files, Stereo-seq GEM, etc.).

The **SpaceTrooper** Bioconductor package (https://bioconductor.org/packages/SpaceTrooper) has already developed robust I/O methods for reading Xenium, CosMx, VisiumHD, and other platform-specific data formats, and its `spatialPerCellQC()` function computes platform-appropriate morphology metrics and stores them in `colData()` of a SpatialExperiment object.

**CellSweeper's approach:**
- CellSweeper operates on **column names in `colData()`**, agnostic to how those columns were computed
- In documentation and vignettes, **recommend that users preprocess their data with SpaceTrooper** (or equivalent) to populate morphology metrics before running CellSweeper
- CellSweeper checks whether specified morphology metric columns exist in `colData()`. If present, they are incorporated into QC. If absent, they are skipped with an informative message
- This keeps CellSweeper's dependency footprint small, avoids duplicating SpaceTrooper's platform-specific parsing logic, and prevents CellSweeper from breaking when platforms change output formats
- SpaceTrooper is listed in `Suggests` (not `Imports`) — it is not a hard dependency

**Example vignette workflow:**
```r
# Step 1: Load and compute morphology metrics (SpaceTrooper)
library(SpaceTrooper)
spe <- readXeniumSPE("/path/to/xenium/output/")
spe <- spatialPerCellQC(spe)

# Step 2: Compute standard QC metrics (scuttle)
library(scuttle)
spe <- addPerCellQCMetrics(spe,
    subsets = list(mito = grep("^MT-", rownames(spe))))

# Step 3: Run CellSweeper
library(CellSweeper)
spe <- runCellSweeper(spe,
    morpho_metrics = c("cell_area", "eccentricity", "compactness"))
```

### Critical Design Decisions

1. **Clustering contamination by QC metrics:** PCA is often driven by library size and related QC metrics. The `clusterCellTypes()` function should support regressing QC covariates out at the gene level BEFORE PCA (not after — regressing from PCs doesn't work well because the damage is baked into the rotation matrix). Regress covariates from the log-normalized expression matrix gene-by-gene, then run PCA on residuals, then cluster. The clustering only needs to be "good enough" for coarse cell types.

2. **Within-cluster spatial neighbors:** When computing k-NN in tissue space within a cluster, some cells may have very few same-cluster spatial neighbors (e.g., rare cell types scattered in the tissue). The `min_neighbors` parameter handles this — if a cell has fewer than `min_neighbors` same-cluster spatial neighbors within a reasonable distance, fall back to global within-cluster QC (not spatial) for that cell.

3. **Spatial dispersion permutation null:** When computing neighborhood homogeneity, larger clusters will have higher homogeneity by chance. The permutation test (randomly permuting all cluster labels, recomputing homogeneity for each cluster) provides a proper null expectation that accounts for this.

4. **Platform generality:** Not all platforms provide the same metrics. Morphological metrics (area, eccentricity) are available for imaging-based platforms (Xenium, CosMx, MERFISH) but may not be available for sequencing-based (Slide-seq). Mito percentage may not be informative for targeted gene panels without mito genes. The API should gracefully handle missing metrics — if a metric column doesn't exist in colData, skip it with a message.

5. **Scalability:** Single-cell resolution datasets can have millions of cells. Use BiocNeighbors (which wraps Annoy/HNSW for approximate NN search) for spatial neighbor computation. The within-cluster NN search is naturally parallelizable across clusters. Consider `BiocParallel` for parallelization.

6. **Transcript density:** Compute as `total_counts / cell_area`. This is a useful derived metric that partially normalizes library size for cell size, helping distinguish truly low-quality cells from small cells with appropriately low counts. Store in colData during the global filtering step.

### What NOT to Implement (Out of Scope for v1)
- Marker gene-based cluster validation (requires a priori knowledge / reference)
- Ambient RNA correction (handled by other tools: SoupX, CellBender, DecontX)
- Doublet detection (handled by: scDblFinder, Scrublet, DoubletFinder)
- Cell segmentation itself (handled by: Cellpose, Baysor, Proseg, 10x built-in)
- Normalization or downstream analysis — CellSweeper is QC only

## Testing Strategy

### Unit Tests
- Test each function independently with small synthetic SPE objects
- Test that colData columns are correctly added
- Test edge cases: clusters with very few cells, cells with no same-cluster spatial neighbors, platforms missing morphological data
- Test that permutation null gives expected results for perfectly random vs perfectly spatially organized clusters

### Integration Tests
- Full pipeline on example Xenium/CosMx datasets
- Verify that flagged cells show expected spatial patterns (clustered near artifacts, at FOV borders, etc.)
- Compare against global MAD-based QC to demonstrate reduced bias across cell types

### Benchmarking (for manuscript, not package tests)
- Compare CellSweeper vs global MAD, ddqc, ctQC, SpotSweeper (naively applied)
- Evaluation metrics:
  - Spatial coherence of clusters before/after QC
  - Downstream clustering quality (silhouette score, ARI with known annotations)
  - Cell-type retention: does CellSweeper retain more legitimate rare cell types?
  - Detection of known artifacts (FOV borders in CosMx, tissue damage regions)
- Datasets: publicly available Xenium, CosMx, Slide-seq, Stereo-seq, VisiumHD

## References

- **SpotSweeper:** Totty M, Hicks SC, Guo B. SpotSweeper: spatially aware quality control for spatial transcriptomics. *Nature Methods* 22, 1520–1530 (2025). https://doi.org/10.1038/s41592-025-02713-3
- **ddqc:** Subramanian A et al. Biology-inspired data-driven quality control for scientific discovery in single-cell transcriptomics. *Genome Biology* 23, 267 (2022).
- **ctQC:** bioRxiv 2024. ctQC improves biological inferences from single cell and spatial transcriptomics data. https://doi.org/10.1101/2024.05.23.594978
- **SpaceTrooper / OSTA imaging QC:** https://lmweber.org/OSTA/pages/img-quality-control.html
- **Library size confounds biology:** Bhuva DD et al. Library size confounds biology in spatial transcriptomics data. *Genome Biology* 25, 99 (2024).
- **OSCA:** Amezquita RA et al. Orchestrating single-cell analysis with Bioconductor. *Nature Methods* 17, 137–145 (2020).
