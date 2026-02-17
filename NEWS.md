# CellSweeper 0.99.0

* Initial Bioconductor submission
* Three-level QC framework for single-cell resolution spatial transcriptomics:
    * Level 1: Global pre-filtering (`globalFilter()`)
    * Level 2: Cluster-level artifact detection (`clusterCellTypes()`, `flagArtifactClusters()`)
    * Level 3: Within-cluster spatial QC (`clusterLocalOutliers()`)
* Convenience wrapper (`runCellSweeper()`)
* Visualization functions (`plotSpatialOutliers()`, `plotClusterQC()`, `plotSpatialDispersion()`, `plotClusterSummary()`)
