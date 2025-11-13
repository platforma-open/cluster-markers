# Overview

Identifies marker genes that distinguish each cell cluster from all other clusters using differential expression analysis. The block performs statistical testing using the Wilcoxon rank-sum test to identify genes significantly enriched in each cluster compared to the rest of the dataset.

Genes are filtered based on configurable thresholds: minimum log2 fold change and maximum adjusted p-value. The strict overlap mode provides additional specificity by requiring that marker genes are expressed in less than 20% of cells in other clusters. The block outputs comprehensive marker gene lists, top markers per cluster for visualization, and differential expression data suitable for functional enrichment analysis.

The identified marker genes can be visualized as dot plots showing expression patterns across clusters, and their spatial distribution can be examined on UMAP or t-SNE embeddings in downstream Cell Browser block. Marker genes are essential for characterizing cell populations, understanding cluster identity, and guiding downstream analyses such as cell type annotation and functional pathway analysis.

The block uses scanpy v1.10.1 for marker gene identification. When using this block in your research, cite the scanpy publication (Wolf et al. 2018) listed below.

The following publication describes the methodology used:

> Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. _Genome Biology_ **19**, 15 (2018). [https://doi.org/10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)
