# Figures

Plots generated automatically by the pipeline scripts.
Not tracked by git — regenerate by running the scripts in order.

## Script 01 — QC and filtering

| Figure | Description |
|--------|-------------|
| `NMU_O_D_vlnplot_QC.png` | Distribution of nFeature_RNA, nCount_RNA, percent.mt and HK_genes across all cells before filtering |
| `NMU_O_P_vlnplot_QC.png` | Same for NMU_O_P |
| `NMU_O_D_scatter_QC.png` | Scatter plots: nCount vs nFeature, nCount vs percent.mt, percent.mt vs HK_genes |
| `NMU_O_P_scatter_QC.png` | Same for NMU_O_P |

## Script 02 — Normalization and dimensionality reduction

| Figure | Description |
|--------|-------------|
| `NMU_O_D_variable_features.png` | Mean-variance plot highlighting the 2,000 most variable genes |
| `NMU_O_P_variable_features.png` | Same for NMU_O_P |
| `NMU_O_D_pca.png` | Cells projected onto PC1 and PC2 |
| `NMU_O_P_pca.png` | Same for NMU_O_P |
| `NMU_O_D_pca_loadings.png` | Top genes contributing to PC1 and PC2 |
| `NMU_O_P_pca_loadings.png` | Same for NMU_O_P |
| `NMU_O_D_elbow.png` | Variance explained per PC — used to select dims = 1:20 |
| `NMU_O_P_elbow.png` | Same for NMU_O_P |
| `NMU_O_D_umap.png` | UMAP embedding after PCA |
| `NMU_O_P_umap.png` | Same for NMU_O_P |
| `NMU_O_D_cellcycle.png` | UMAP colored by predicted cell cycle phase (G1, S, G2M) |
| `NMU_O_P_cellcycle.png` | Same for NMU_O_P |
| `NMU_O_D_doublets.png` | UMAP colored by DoubletFinder classification (Singlet/Doublet) |
| `NMU_O_P_doublets.png` | Same for NMU_O_P |

## Script 03 — Clustering and annotation

| Figure | Description |
|--------|-------------|
| `NMU_O_D_clustree.png` | Cluster stability across resolution sweep (0 to 1). Resolution 0.1 selected |
| `NMU_O_P_clustree.png` | Same for NMU_O_P. Resolution 0.2 selected |
| `NMU_O_D_dotplot_markers.png` | Top 5 markers per cluster: dot size = detection rate, color = expression level |
| `NMU_O_P_dotplot_markers.png` | Same for NMU_O_P |
| `NMU_O_D_featureplot_canonical.png` | Expression of Krt5 (basal) and Upk3a (luminal) on UMAP |
| `NMU_O_P_featureplot_canonical.png` | Expression of Mki67 (proliferation), Krt14 (basal) and Psca (intermediate) on UMAP |
| `NMU_O_D_umap_annotated.png` | Final UMAP: Intermediate, Basal, Luminal |
| `NMU_O_P_umap_annotated.png` | Final UMAP: Basal, Intermediate high, Intermediate low, Basal G2M |

## Script 04 — Integration and downstream analysis

| Figure | Description |
|--------|-------------|
| `merged_pca_by_sample.png` | PCA of merged object colored by sample — shows batch effect before correction |
| `integrated_umap_by_sample.png` | Harmony-corrected UMAP colored by sample — cells should be interleaved |
| `integrated_umap_by_cluster.png` | Harmony-corrected UMAP colored by annotated cell population |
| `population_composition.png` | Stacked bar chart: proportion of each cell type in NMU_O_D vs NMU_O_P |
