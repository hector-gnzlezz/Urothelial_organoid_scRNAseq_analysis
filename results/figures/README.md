# Figures

Plots generated automatically by the pipeline scripts.

## Script 01 — QC and filtering

| Figure | Description |
|--------|-------------|
| `01/vln_plot_01.png` | Distribution of nFeature_RNA, nCount_RNA, percent.mt and HK_genes across all cells before filtering (NMU_O_D) |
| `01/vln_plot_02.png` | Same for NMU_O_P |
| `01/quality_HK_mt_NMU_O_D.png` | Scatter plots: nCount vs nFeature, nCount vs percent.mt, percent.mt vs HK_genes (NMU_O_D) |
| `01/quality_HK_mt_NMU_O_P.png` | Same for NMU_O_P |

## Script 02 — Normalization and dimensionality reduction

| Figure | Description |
|--------|-------------|
| `02/VariableFeaturePlot(NMU_O_D).png` | Mean-variance plot highlighting the 2,000 most variable genes (NMU_O_D) |
| `02/VariableFeaturePlot(NMU_O_P).png` | Same for NMU_O_P |
| `02/PCA.png` | Cells projected onto PC1 and PC2 |
| `02/Top_genes_D.png` | Top genes contributing to PC1 and PC2 (NMU_O_D) |
| `02/Top_genes_P.png` | Same for NMU_O_P |
| `02/JackStraw_D.png` | Statistical significance of each PC (NMU_O_D) |
| `02/JackStraw_P.png` | Same for NMU_O_P |
| `02/ElbowPlot_D.png` | Variance explained per PC — used to select dims = 1:20 (NMU_O_D) |
| `02/ElbowPlot_P.png` | Same for NMU_O_P |
| `02/Heatmap_D.png` | Top genes driving each of the first 15 PCs (NMU_O_D) |
| `02/Heatmap_P.png` | Same for NMU_O_P |
| `02/UMAP_D.png` | UMAP embedding after PCA (NMU_O_D) |
| `02/UMAP_P.png` | Same for NMU_O_P |
| `02/UMAP_D_cell_cycle.png` | UMAP colored by predicted cell cycle phase: G1, S, G2M (NMU_O_D) |
| `02/UMAP_P_cell_cycle.png` | Same for NMU_O_P |
| `02/UMAP_D_DoubletFinder.png` | UMAP colored by DoubletFinder classification: Singlet/Doublet (NMU_O_D) |
| `02/UMAP_P_DoubletFinder.png` | Same for NMU_O_P |
| `02/BCmetric_DoubletFinder_D.png` | pK parameter sweep for DoubletFinder (NMU_O_D) |
| `02/BCmetric_DoubletFinder_P.png` | Same for NMU_O_P |

## Script 03 — Clustering and annotation

| Figure | Description |
|--------|-------------|
| `03/Clustree(NMU_O_D).png` | Cluster stability across resolution sweep 0–1. Resolution 0.1 selected (NMU_O_D) |
| `03/Clustree(NMU_O_P).png` | Same for NMU_O_P. Resolution 0.2 selected |
| `03/Dimplot_D_clustered.png` | UMAP colored by cluster before annotation (NMU_O_D) |
| `03/Dimplot_P_clustered.png` | Same for NMU_O_P |
| `03/DotPlot_TOP_GeneExpression_D.png` | Top 5 markers per cluster: dot size = detection rate, color = expression level (NMU_O_D) |
| `03/DotPlot_TOP_GeneExpression_P.png` | Same for NMU_O_P |
| `03/Krt5_Upk3a_Expression_D.png` | Expression of Krt5 (basal) and Upk3a (luminal) on UMAP (NMU_O_D) |
| `03/Mki67_Krt14_Psca_expression_P.png` | Expression of Mki67 (proliferation), Krt14 (basal) and Psca (intermediate) on UMAP (NMU_O_P) |
| `03/Cell_Type_D.png` | Final annotated UMAP: Intermediate, Basal, Luminal (NMU_O_D) |
| `03/Cell_Type_P.png` | Final annotated UMAP: Basal, Intermediate high, Intermediate low, Basal G2M (NMU_O_P) |

## Script 04 — Integration and downstream analysis

| Figure | Description |
|--------|-------------|
| `04/PCA_Merged.png` | PCA of merged object colored by sample — shows batch effect before correction |
| `04/DimPlot_Color_by_Sample.png` | Harmony-corrected UMAP colored by sample — cells should be interleaved |
| `04/DimPlot_Color_by_Custer.png` | Harmony-corrected UMAP colored by annotated cell population |
| `04/ClusTree_Merged.png` | Cluster stability in the integrated object across resolution sweep |
| `04/Annotated_Cell_Type.png` | Final annotated UMAP of the integrated dataset |
| `04/Population_Plot.png` | Stacked bar chart: proportion of each cell type in NMU_O_D vs NMU_O_P |
