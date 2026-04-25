# Project: Urothelial organoid scRNA-seq analysis
# Script: 04 - Sample integration and downstream analysis
# Description: Merge samples, correct batch effect with Harmony,
#              integrated clustering and gene signature export for ssGSEA/CIBERSORTx
# Input: Annotated Seurat objects from script 03 (NMU_O_D, NMU_O_P)

library(harmony)    # Batch correction via iterative PCA alignment
library(dplyr)      # Pipe operator and data manipulation
library(ggplot2)    # Visualization
library(nichenetr)  # Mouse-to-human gene symbol conversion


# 1. MERGE SAMPLES

# Before integration we need to merge both Seurat objects into one.
# add.cell.ids adds a prefix to each cell barcode ("D_" or "P_")
# to avoid barcode collisions between the two datasets
# (10X barcodes are not unique across samples)
NMU_O_merged <- merge(x = NMU_O_D, y = NMU_O_P, add.cell.ids = c("D", "P"))
NMU_O_merged

# Seurat v5 stores layers separately per sample after merging.
# JoinLayers collapses them into a single counts layer,
# which is required before running normalization on the merged object
NMU_O_merged <- JoinLayers(NMU_O_merged)


# 2. PREPROCESSING OF MERGED OBJECT

# Repeat the full preprocessing pipeline on the merged object.
# This is necessary because normalization and PCA must be computed
# on the combined dataset to place all cells in the same feature space
# before batch correction can be applied.
NMU_O_merged <- NormalizeData(NMU_O_merged, normalization.method = "LogNormalize", scale.factor = 10000)
NMU_O_merged <- FindVariableFeatures(NMU_O_merged, selection.method = "vst", nfeatures = 2000)
NMU_O_merged <- ScaleData(NMU_O_merged)
NMU_O_merged <- RunPCA(NMU_O_merged, features = VariableFeatures(object = NMU_O_merged))

# Visualize merged PCA colored by sample of origin.
# If samples separate strongly along PC1/PC2, batch effect is present
# and Harmony correction is needed before clustering.
DimPlot(NMU_O_merged, reduction = "pca", group.by = "orig.ident")


# 3. BATCH CORRECTION WITH HARMONY

# Harmony iteratively adjusts PCA embeddings to remove variation
# attributable to the batch variable (here: sample of origin, "orig.ident").
# It works in PC space — it does not modify the count matrix.
# The output is a corrected "harmony" dimensional reduction that can be used
# instead of "pca" for UMAP and clustering.
# This is critical here because NMU_O_D and NMU_O_P come from different
# culture conditions — without correction, cells would cluster by condition
# rather than by cell type.
NMU_O_integrated <- RunHarmony(NMU_O_merged, "orig.ident")

# Run UMAP on the Harmony-corrected embeddings (not PCA)
NMU_O_integrated <- NMU_O_integrated %>%
  RunUMAP(reduction = "harmony", dims = 1:20)


# 4. CLUSTERING ON INTEGRATED DATA

# FindNeighbors and FindClusters use the Harmony reduction instead of PCA,
# ensuring that clustering reflects cell type rather than batch
NMU_O_integrated <- NMU_O_integrated %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

# Select resolution using clustree — same logic as script 03
clustree(NMU_O_integrated)
Idents(NMU_O_integrated) <- NMU_O_integrated$RNA_snn_res.0.1

# Two key visualizations after integration:
# 1. Color by sample — cells should be interleaved, not separated by condition.
#    If separation persists, Harmony correction was insufficient.
DimPlot(NMU_O_integrated, reduction = "umap", group.by = "orig.ident")

# 2. Color by cluster — biological populations should be shared across samples
DimPlot(NMU_O_integrated, reduction = "umap")


# 5. POPULATION COMPOSITION ANALYSIS

# Stacked bar chart showing the proportion of each cell population
# in the differentiated vs proliferating condition.
# This is the key biological result: how do population proportions
# shift between conditions?
# Expected: proliferating condition enriched in Basal G2M,
#           differentiated condition enriched in Luminal cells
plot_data <- as.data.frame(table(Idents(NMU_O_integrated), NMU_O_integrated$orig.ident))
colnames(plot_data) <- c("Cluster", "Sample", "Count")

ggplot(plot_data, aes(x = Sample, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(fill = "Cell population") +
  theme_minimal() +
  theme(
    text         = element_text(size = 18),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.title   = element_blank(),
    axis.text.y  = element_blank()
  )


# 6. ANNOTATION AND MARKER EXPORT

# Annotate integrated clusters using the same biological logic as script 03.
# Always verify cluster order with levels(NMU_O_integrated) before running.
new.cluster.ids <- c("Intermediate", "Basal", "Luminal", "Basal G2M")
names(new.cluster.ids) <- levels(NMU_O_integrated)
NMU_O_integrated <- RenameIdents(NMU_O_integrated, new.cluster.ids)
DimPlot(NMU_O_integrated, reduction = "umap", label = TRUE)

# Find markers for all populations in the integrated dataset
# These will be used as gene signatures for downstream bulk RNA-seq analysis
NMU_O_integrated.markers <- FindAllMarkers(
  NMU_O_integrated,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

# Export full marker table for reference
write.table(
  NMU_O_integrated.markers,
  "./results/NMU_O_integrated_markers.txt",
  sep = '\t'
)

# Save full R session to allow resuming analysis without re-running all scripts
save.image(file = "./results/scRNAseq.RData")


# 7. GENE SIGNATURE EXPORT FOR ssGSEA

# Strategy: use scRNA-seq cell population markers as gene signatures
# to score bulk RNA-seq samples via single-sample GSEA (ssGSEA).
# This "projects" the single-cell findings onto bulk data,
# estimating the relative abundance of each cell population per bulk sample.
# Tool used: GenePattern ssGSEA module (https://cloud.genepattern.org)

# Filter markers to statistically significant genes only
significant_markers <- NMU_O_integrated.markers[NMU_O_integrated.markers$p_val_adj < 0.05, ]

# Build a named list of gene vectors, one per cell population
gene_sets <- lapply(unique(significant_markers$cluster), function(cluster_id) {
  cluster_markers <- significant_markers[significant_markers$cluster == cluster_id, ]
  return(cluster_markers$gene)
})
names(gene_sets) <- unique(NMU_O_integrated.markers$cluster)

# Convert mouse gene symbols to human orthologs.
# ssGSEA will be run on human bulk RNA-seq data, so gene symbols must match.
# nichenetr uses a curated ortholog table for this conversion.
Basal_G2M_GS  <- as.data.frame(convert_mouse_to_human_symbols(gene_sets[["Basal G2M"]],  version = 1))
Basal_GS       <- as.data.frame(convert_mouse_to_human_symbols(gene_sets[["Basal"]],       version = 1))
Intermediate_GS <- as.data.frame(convert_mouse_to_human_symbols(gene_sets[["Intermediate"]], version = 1))
Luminal_GS     <- as.data.frame(convert_mouse_to_human_symbols(gene_sets[["Luminal"]],     version = 1))

# Rebuild gene set list with human symbols
gene_sets <- list(
  Basal_G2M_GS[,1],
  Basal_GS[,1],
  Intermediate_GS[,1],
  Luminal_GS[,1]
)
names(gene_sets) <- unique(NMU_O_integrated.markers$cluster)

# GMT format: each line = one gene set
# Format: <set_name> \t <description> \t <gene1> \t <gene2> \t ...
# NA values (genes without human orthologs) are removed before export
save_as_gmt <- function(gene_sets, file_name) {
  file_conn <- file(file_name, open = "w")
  for (set_name in names(gene_sets)) {
    cleaned_genes <- na.omit(gene_sets[[set_name]])
    if (length(cleaned_genes) > 0) {
      line <- paste(set_name, "description", paste(cleaned_genes, collapse = "\t"), sep = "\t")
      writeLines(line, file_conn)
    }
  }
  close(file_conn)
}

save_as_gmt(gene_sets, "./results/NMU_O_gene_sets.gmt")


# 8. REFERENCE MATRIX EXPORT FOR CIBERSORTx DECONVOLUTION

# Strategy: use a small subset of annotated scRNA-seq cells as a reference
# to deconvolute bulk RNA-seq samples via CIBERSORTx.
# CIBERSORTx estimates the fraction of each cell population in each bulk sample
# using the single-cell expression profiles as a reference signature matrix.
# Tool: https://cibersortx.stanford.edu/

# Sample 10 cells per population (minimum required by CIBERSORTx)
# set.seed ensures the same cells are selected if the script is re-run
set.seed(123)
idents <- unique(Idents(NMU_O_integrated))
subset_cells <- unlist(lapply(idents, function(ident) {
  sample(
    Cells(NMU_O_integrated)[Idents(NMU_O_integrated) == ident],
    size    = min(10, sum(Idents(NMU_O_integrated) == ident)),
    replace = FALSE
  )
}))

NMU_O_subset <- subset(NMU_O_integrated, cells = subset_cells)

# Extract raw count matrix (CIBERSORTx expects counts, not normalized values)
NMU_O_expr_matrix <- as.data.frame(GetAssayData(NMU_O_subset, slot = "counts"))

# Label columns with cell population names (required by CIBERSORTx format)
colnames(NMU_O_expr_matrix) <- Idents(NMU_O_subset)

# Convert mouse gene symbols to human orthologs (bulk data uses human symbols)
# Remove genes without orthologs (NA) and duplicated conversions
human_converted <- as.data.frame(
  convert_mouse_to_human_symbols(rownames(NMU_O_expr_matrix), version = 1)
)
NMU_O_expr_matrix$human_converted <- human_converted[, 1]
NMU_O_expr_matrix <- na.omit(NMU_O_expr_matrix)
NMU_O_expr_matrix <- NMU_O_expr_matrix[!duplicated(NMU_O_expr_matrix$human_converted), ]
rownames(NMU_O_expr_matrix) <- NMU_O_expr_matrix$human_converted
NMU_O_expr_matrix$human_converted <- NULL

# Export as tab-separated matrix with row names (gene symbols)
write.table(
  NMU_O_expr_matrix,
  file      = "./results/NMU_O_expr_matrix.tsv",
  sep       = '\t',
  col.names = NA
)
