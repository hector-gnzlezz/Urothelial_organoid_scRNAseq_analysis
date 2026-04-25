# Project: Urothelial organoid scRNA-seq analysis
# Script: 03 - Clustering and cell type annotation
# Description: Resolution sweep with clustree, marker identification,
#              canonical marker visualization and population annotation
# Input: Seurat objects after normalization and dim reduction (scripts 01-02)

library(clustree) # For resolution sweep visualization


# 1. GRAPH-BASED CLUSTERING

# Seurat uses a graph-based clustering approach (Louvain algorithm):
# 1. FindNeighbors builds a KNN graph in PC space — each cell is connected
#    to its k nearest neighbors based on euclidean distance in PCA coordinates
# 2. FindClusters partitions this graph by optimizing modularity
#    (how densely connected cells are within clusters vs between clusters)
# This is the "Louvain" unsupervised clustering method mentioned in the TFM context
# as a "black box" — it groups cells but does not explain why

# dims = 1:20 matches the number of significant PCs selected in script 02
NMU_O_D <- FindNeighbors(NMU_O_D, dims = 1:20)
NMU_O_P <- FindNeighbors(NMU_O_P, dims = 1:20)

# Resolution controls cluster granularity:
# low resolution (0.1) = few large clusters
# high resolution (1.0) = many small clusters
# We sweep across a range to find the most biologically meaningful resolution
NMU_O_D <- FindClusters(NMU_O_D, resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
NMU_O_P <- FindClusters(NMU_O_P, resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))


# 2. RESOLUTION SELECTION WITH CLUSTREE

# Clustree visualizes how clusters split as resolution increases.
# A stable cluster (consistent across resolutions) is more likely to represent
# a real biological population rather than an artifact of over-clustering.
# We look for the lowest resolution where biologically meaningful populations
# are separated — adding more resolution should not split stable clusters.
clustree(NMU_O_D)
clustree(NMU_O_P)

# Based on clustree inspection:
# NMU_O_D: resolution 0.1 gives 3 stable clusters matching known urothelial populations
# NMU_O_P: resolution 0.2 gives 4 stable clusters, capturing the extra proliferating population
Idents(NMU_O_D) <- NMU_O_D$RNA_snn_res.0.1
Idents(NMU_O_P) <- NMU_O_P$RNA_snn_res.0.2

# Visualize selected clustering on UMAP
DimPlot(NMU_O_D, reduction = "umap")
DimPlot(NMU_O_P, reduction = "umap")


# 3. DIFFERENTIAL MARKER IDENTIFICATION

# FindAllMarkers runs a Wilcoxon rank-sum test for each cluster vs all other cells,
# identifying genes that are significantly overexpressed in each cluster.
# These markers are used to biologically annotate each cluster.
# Filters applied:
# - only.pos = TRUE: keep only upregulated markers (positive markers are more interpretable)
# - min.pct = 0.25: gene must be expressed in at least 25% of cells in the cluster
# - logfc.threshold = 0.25: minimum log2 fold-change (removes very subtle differences)
NMU_O_D_0.1.markers <- FindAllMarkers(
  NMU_O_D,
  only.pos       = TRUE,
  min.pct        = 0.25,
  logfc.threshold = 0.25
)

NMU_O_P_0.2.markers <- FindAllMarkers(
  NMU_O_P,
  only.pos       = TRUE,
  min.pct        = 0.25,
  logfc.threshold = 0.25
)

# Extract top 5 markers per cluster ranked by average log2 fold-change
# These are the genes most specifically expressed in each cluster
NMU_O_D_0.1.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5_D

NMU_O_P_0.2.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5_P


# 4. MARKER VISUALIZATION

# DotPlot shows two dimensions simultaneously for each gene per cluster:
# - dot size = percentage of cells expressing the gene (detection rate)
# - dot color = average expression level (blue = low, red = high)
# This allows quick identification of cluster-specific marker genes
DotPlot(NMU_O_D, features = top5_D$gene) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradientn(colors = c("blue", "grey", "red"))

DotPlot(NMU_O_P, features = top5_P$gene) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradientn(colors = c("blue", "grey", "red"))

# FeaturePlot overlays gene expression on UMAP to validate canonical marker identity.
# This is the biological ground truth that confirms our cluster annotations:
# NMU_O_D canonical markers:
# - Krt5: basal cell marker (cytokeratin expressed in progenitor/basal layer)
# - Upk3a: luminal/umbrella cell marker (uroplakin, terminal differentiation)
FeaturePlot(NMU_O_D, features = c("Krt5", "Upk3a"), cols = c("grey", "red"), reduction = "umap")

# NMU_O_P canonical markers:
# - Mki67: proliferation marker (Ki67, expressed in actively dividing cells)
# - Krt14: basal progenitor marker (cytokeratin 14)
# - Psca: intermediate/luminal progenitor marker (prostate stem cell antigen)
FeaturePlot(NMU_O_P, features = c("Mki67", "Krt14", "Psca"), cols = c("grey", "red"), reduction = "umap")


# 5. CELL TYPE ANNOTATION
# The annotation step at the end is the most biologically critical part of the whole pipeline and also the most manual. The cluster numbers Seurat assigns are arbitrary, cluster 0 is just the largest cluster, not necessarily "Basal". The annotation only works if the order in new.cluster.ids exactly matches levels(NMU_O_D). Always run levels(NMU_O_D) first and double-check the order before renaming. 

# Based on the combination of:
# 1. Data-driven markers from FindAllMarkers
# 2. Canonical marker expression from FeaturePlot
# 3. Known urothelial biology (basal → intermediate → luminal differentiation axis)
# We assign biological identities to each cluster

# NMU_O_D (differentiated condition): 3 populations
# Cluster order must match levels(NMU_O_D) — verify with levels() before running
new.cluster.ids <- c("Intermediate", "Basal", "Luminal")
names(new.cluster.ids) <- levels(NMU_O_D)
NMU_O_D <- RenameIdents(NMU_O_D, new.cluster.ids)
DimPlot(NMU_O_D, reduction = "umap", label = TRUE)

# NMU_O_P (proliferating condition): 4 populations
# Extra population vs differentiated: Basal G2M (actively cycling basal cells)
# This is consistent with the proliferating culture condition
new.cluster.ids <- c("Basal", "Intermediate high", "Intermediate low", "Basal G2M")
names(new.cluster.ids) <- levels(NMU_O_P)
NMU_O_P <- RenameIdents(NMU_O_P, new.cluster.ids)
DimPlot(NMU_O_P, reduction = "umap", label = TRUE)
