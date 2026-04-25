# Project: Urothelial organoid scRNA-seq analysis
# Script: 02 - Normalization and dimensionality reduction
# Description: Normalize counts, identify variable features, run PCA/UMAP,
#              score cell cycle phases and detect doublets
# Input: Filtered Seurat objects from script 01 (NMU_O_D, NMU_O_P)

library(nichenetr)   # For human-to-mouse gene symbol conversion
library(DoubletFinder) # For computational doublet detection


# 1. NORMALIZATION

# Log normalization: divides each cell's counts by the total counts for that cell,
# multiplies by a scale factor (10,000 = counts per 10k), then log-transforms.
# This corrects for differences in sequencing depth between cells:
# a gene expressed in 10 out of 100 total counts is comparable to
# the same gene expressed in 1000 out of 10,000 total counts
NMU_O_D <- NormalizeData(NMU_O_D, normalization.method = "LogNormalize", scale.factor = 10000)
NMU_O_P <- NormalizeData(NMU_O_P, normalization.method = "LogNormalize", scale.factor = 10000)


# 2. IDENTIFICATION OF HIGHLY VARIABLE FEATURES

# Not all genes are informative for distinguishing cell types.
# FindVariableFeatures selects genes with high cell-to-cell variation
# (highly expressed in some cells, lowly expressed in others).
# These 2000 genes will drive the downstream PCA and clustering.
# Method "vst" (variance stabilizing transformation) models the mean-variance
# relationship and selects genes that deviate from the expected trend.
NMU_O_D <- FindVariableFeatures(NMU_O_D, selection.method = "vst", nfeatures = 2000)
NMU_O_P <- FindVariableFeatures(NMU_O_P, selection.method = "vst", nfeatures = 2000)

# Inspect the top 10 most variable genes as a sanity check:
# in urothelial tissue we expect to see differentiation markers here
top10_D <- head(VariableFeatures(NMU_O_D), 10)
top10_P <- head(VariableFeatures(NMU_O_P), 10)

# Plot all variable features; label the top 10
# Genes in the upper right are highly variable and will be most informative
plot1 <- VariableFeaturePlot(NMU_O_D)
plot2 <- LabelPoints(plot = plot1, points = top10_D, repel = TRUE)
plot2

plot1 <- VariableFeaturePlot(NMU_O_P)
plot2 <- LabelPoints(plot = plot1, points = top10_P, repel = TRUE)
plot2


# 3. SCALING

# ScaleData z-scores each gene across all cells (mean = 0, variance = 1).
# This prevents highly expressed genes from dominating PCA simply because
# of their absolute expression level, not because of their biological variation.
# Only applied to variable features by default.
NMU_O_D <- ScaleData(NMU_O_D)
NMU_O_P <- ScaleData(NMU_O_P)


# 4. LINEAR DIMENSIONALITY REDUCTION (PCA)

# PCA reduces the 2000 variable features to a smaller set of principal components (PCs)
# that capture the main axes of variation in the data.
# Each PC is a linear combination of genes; the first PCs capture the most variance.
# This is the same PCA concept from your ML lectures, applied to gene expression.
NMU_O_D <- RunPCA(NMU_O_D, features = VariableFeatures(object = NMU_O_D))
NMU_O_P <- RunPCA(NMU_O_P, features = VariableFeatures(object = NMU_O_P))

# Visualize cells in PCA space — distinct clusters suggest real biological structure
DimPlot(NMU_O_D, reduction = "pca")
DimPlot(NMU_O_P, reduction = "pca")

# VizDimLoadings shows which genes contribute most to each PC.
# High loadings on PC1/PC2 indicate genes that drive the main biological variation.
# In urothelial tissue, expect differentiation markers (Upk3a, Krt5) to appear here.
VizDimLoadings(NMU_O_D, dims = 1:2, reduction = "pca")
VizDimLoadings(NMU_O_P, dims = 1:2, reduction = "pca")

# Heatmaps show the top genes driving each PC across 500 cells.
# A clear gradient from high to low expression confirms the PC captures
# real biological variation rather than noise.
DimHeatmap(NMU_O_D, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(NMU_O_P, dims = 1:15, cells = 500, balanced = TRUE)


# 5. SELECT NUMBER OF SIGNIFICANT PCs

# We need to decide how many PCs to use for downstream clustering and UMAP.
# Using too few loses biological signal; using too many adds noise.
# Two complementary approaches:

# JackStraw: permutation test that identifies statistically significant PCs.
# Computationally expensive — comment out if running on a laptop.
# PCs with a significant p-value (low p, high score) should be retained.
NMU_O_D <- JackStraw(NMU_O_D, num.replicate = 100, dims = 30)
NMU_O_D <- ScoreJackStraw(NMU_O_D, dims = 1:30)
JackStrawPlot(NMU_O_D, dims = 1:30)

# ElbowPlot: faster heuristic. Look for the "elbow" where variance explained
# stops dropping sharply — PCs after the elbow add little information.
# Here we use dims = 1:20 based on the elbow observed in the plot.
ElbowPlot(NMU_O_D, ndims = 30)

NMU_O_P <- JackStraw(NMU_O_P, num.replicate = 100, dims = 30)
NMU_O_P <- ScoreJackStraw(NMU_O_P, dims = 1:30)
JackStrawPlot(NMU_O_P, dims = 1:30)
ElbowPlot(NMU_O_P, ndims = 30)


# 6. NON-LINEAR DIMENSIONALITY REDUCTION (UMAP)

# UMAP projects the high-dimensional PC space into 2D for visualization.
# Unlike PCA, UMAP preserves local neighborhood structure:
# cells that are similar to each other in PC space will be close in UMAP space.
# We use dims = 1:20 based on the elbow plot above.
# Note: UMAP is for visualization only — clustering is done in PC space, not UMAP.
NMU_O_D <- RunUMAP(NMU_O_D, dims = 1:20)
DimPlot(NMU_O_D, reduction = "umap")

NMU_O_P <- RunUMAP(NMU_O_P, dims = 1:20)
DimPlot(NMU_O_P, reduction = "umap")


# 7. CELL CYCLE SCORING

# Cell cycle phase can be a major source of transcriptional variation,
# potentially masking biologically relevant differences between cell types.
# Scoring each cell allows us to visualize and account for cell cycle effects.

# Seurat provides human cell cycle marker genes built-in (cc.genes).
# Since our data is mouse, we convert human gene symbols to mouse orthologs
# using nichenetr's conversion function.
s.genes   <- cc.genes$s.genes    # S phase markers
g2m.genes <- cc.genes$g2m.genes  # G2/M phase markers

s.genes   <- convert_human_to_mouse_symbols(s.genes,   version = 1)
g2m.genes <- convert_human_to_mouse_symbols(g2m.genes, version = 1)

# Assign each cell a cell cycle score and predicted phase (G1, S, or G2M)
NMU_O_D <- CellCycleScoring(NMU_O_D, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(NMU_O_D@meta.data) # Check that S.Score, G2M.Score and Phase columns were added
DimPlot(NMU_O_D, group.by = "Phase", reduction = "umap")
# If cells cluster by phase rather than cell type, consider regressing out
# cell cycle scores in ScaleData (vars.to.regress = c("S.Score", "G2M.Score"))

NMU_O_P <- CellCycleScoring(NMU_O_P, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(NMU_O_P@meta.data)
DimPlot(NMU_O_P, group.by = "Phase", reduction = "umap")


# 8. DOUBLET DETECTION (DoubletFinder)

# Doublets are droplets that captured two cells instead of one.
# They appear as artificial "intermediate" cell types and can confound clustering.
# DoubletFinder detects them by:
# 1. Generating artificial doublets by averaging pairs of real cells
# 2. Embedding real + artificial cells together in PC space
# 3. Scoring each real cell by how many of its nearest neighbors are artificial
#    (pANN score — proportion of Artificial Nearest Neighbors)
# 4. Classifying the top N cells by pANN score as doublets,
#    where N = expected number of doublets based on capture rate

# --- NMU_O_D ---
# 10X captures ~5% doublets per 10,000 cells loaded, adjust if your loading differs
nExp <- round(ncol(NMU_O_D) * 0.05)

# pK is the PC neighborhood size parameter. We optimize it empirically
# by sweeping across values and choosing the one with the highest
# BCmetric (bimodality coefficient — doublets should form a distinct population)
sweep.res  <- paramSweep(NMU_O_D, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn      <- find.pK(sweep.stats)
best_pK    <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Run DoubletFinder with optimized pK
# pN = 0.25: ratio of artificial doublets to real cells (standard default)
NMU_O_D <- doubletFinder(NMU_O_D, PCs = 1:20, pN = 0.25, pK = best_pK, nExp = nExp, sct = FALSE)

# The classification column name includes pN, pK and nExp values 
# update this string if your parameters differ from the example below
# column name will be different every time DoubletFinder is run, because it encodes your specific pK and nExp values. Its necessary to update it after each run
NMU_O_D$DF_classification <- NMU_O_D$DF.classifications_0.25_0.27_239

# Visualize doublets on UMAP — they should not cluster cleanly with any population
DimPlot(NMU_O_D, group.by = "DF_classification")
table(NMU_O_D$DF_classification) # Check absolute numbers of singlets vs doublets


# --- NMU_O_P ---
nExp <- round(ncol(NMU_O_P) * 0.05)

sweep.res   <- paramSweep(NMU_O_P, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn       <- find.pK(sweep.stats)
best_pK     <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

NMU_O_P <- doubletFinder(NMU_O_P, PCs = 1:20, pN = 0.25, pK = best_pK, nExp = nExp, sct = FALSE)

NMU_O_P$DF_classification <- NMU_O_P$DF.classifications_0.25_0.01_385

DimPlot(NMU_O_P, group.by = "DF_classification")
table(NMU_O_P$DF_classification)
