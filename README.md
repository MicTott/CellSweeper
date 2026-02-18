
**Cluster-aware spatial quality control for single-cell resolution spatial transcriptomics.**

CellSweeper provides a QC framework for single-cell resolution spatial transcriptomics platforms such as Xenium, CosMx, MERFISH, VisiumHD, and others. It extends the [SpotSweeper](https://bioconductor.org/packages/SpotSweeper) local outlier framework in control for confounding biology (Totty, Hicks & Guo, *Nature Methods* 2025), while also incorporating segmentation morphology metrics from [SpaceTrooper](https://github.com/drighelli/SpaceTrooper). 

## The Three-Level QC Framework

| Level | Function(s) | What it does |
|-------|-------------|--------------|
| **1. Global pre-filtering** | `globalFilter()` | Remove obvious low quality observations — near-zero counts, impossible segmentations, etc |
| **2. Cluster-level QC** | `clusterCellTypes()` + `flagArtifactClusters()` | Cluster cells in gene expression space, then flag entire artifact clusters via multivariate pseudobulk outlier detection and spatial dispersion analysis |
| **3. Within-cluster local outlier detection** | `clusterLocalOutliers()` | Local outlier detection restricted to same-cluster neighbors, applied to both transcriptomic and morphological metrics |

## Installation

```r
## Install from Bioconductor (coming soon)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("CellSweeper")

# Or install the development version from GitHub
BiocManager::install("MicTott/CellSweeper")
```

## Integration with SpaceTrooper

For imaging-based platforms (Xenium, CosMx, MERFISH), we recommend preprocessing with [SpaceTrooper](https://bioconductor.org/packages/SpaceTrooper) to compute morphology-based QC metrics and cell boundary polygons before running CellSweeper:

```r
library(SpaceTrooper)
spe <- readXeniumSPE("/path/to/xenium/output/")
spe <- spatialPerCellQC(spe)    # adds Area_um, log2AspectRatio, etc.
spe <- addPolygonsToSPE(spe)    # adds sf polygons for visualization

library(CellSweeper)
spe <- runCellSweeper(spe,
    morpho_metrics = c("Area_um", "log2AspectRatio", "log2CountArea"))
```

## Citation

If you use CellSweeper in your research, please cite:

> Totty M, Hicks SC, Guo B. CellSweeper: spatially-aware quality control for single cell spatial transcriptomics. *In preparation.*

If you use SpaceTrooper for preprocessing, please also cite:

>  Add SpaceTrooper citation here:

CellSweeper is a direct extension of SpotSweeper:

> Totty M, Hicks SC, Guo B. SpotSweeper: spatially aware quality control for spatial transcriptomics. *Nature Methods* 22, 1520–1530 (2025). [doi:10.1038/s41592-025-02713-3](https://doi.org/10.1038/s41592-025-02713-3)

## Report Bugs and Suggest Features.

Please open an [issue](https://github.com/MicTott/CellSweeper/issues) to report bugs or suggest features.
