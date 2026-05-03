# Project: Urothelial organoid scRNA-seq analysis
# Script: 01 - QC and filtering
# Description: Load 10X data, compute QC metrics and filter low-quality cells

# Set working directory to the root of the repository
# All paths below are relative to this root
library(here)

# --- Libraries ---
library(Seurat)    # Core scRNA-seq analysis framework
library(dplyr)     # Data manipulation
library(patchwork) # Combine multiple ggplots into one figure
library(ggplot2)   # Visualization


# 1. LOAD DATA

# Read10X loads the three files output by Cell Ranger (barcodes, features, matrix)
# and returns a sparse count matrix (genes x cells)
NMU_O_D.data <- Read10X(data.dir = here("data/raw/NMU_O_D/"))

# CreateSeuratObject wraps the count matrix into a Seurat object
# min.cells = 3: keep only genes detected in at least 3 cells (removes noise)
# min.features = 200: keep only cells with at least 200 detected genes (removes empty droplets)
NMU_O_D <- CreateSeuratObject(
  counts       = NMU_O_D.data,
  project      = "NMU_O_D",
  min.cells    = 3,
  min.features = 200
)
NMU_O_D # Print summary: number of genes and cells retained


# --- Why sparse matrices matter ---
# 10X data is mostly zeros (most genes are not expressed in most cells)
# Storing as a dense matrix would use ~10x more memory
# This comparison shows the practical difference:
dense.size  <- object.size(as.matrix(NMU_O_D.data)) # Convert to dense and measure
sparse.size <- object.size(NMU_O_D.data)            # Measure sparse matrix
dense.size   # Will be much larger
sparse.size  # Efficient sparse format

# Check expression of known urothelial marker genes across first 30 cells
# Upk3a = luminal marker, Krt5 = basal marker, Psca = intermediate marker
# Most values will be 0, illustrating data sparsity
NMU_O_D.data[c("Upk3a", "Krt5", "Psca"), 1:30]


# Repeat for the proliferating condition (NMU_O_P)
NMU_O_P.data <- Read10X(data.dir = here("data/raw/NMU_O_P/"))
NMU_O_P <- CreateSeuratObject(
  counts       = NMU_O_P.data,
  project      = "NMU_O_P",
  min.cells    = 3,
  min.features = 200
)
NMU_O_P


# 2. COMPUTE QC METRICS

# --- Mitochondrial gene percentage ---
# High mitochondrial gene expression (>20%) indicates damaged or dying cells:
# when a cell lyses, cytoplasmic mRNA leaks out but mitochondrial mRNA is retained,
# so the remaining reads are disproportionately mitochondrial
# Mouse mitochondrial genes are prefixed with "mt-" (lowercase)
NMU_O_D[["percent.mt"]] <- PercentageFeatureSet(NMU_O_D, pattern = "^mt-")
NMU_O_P[["percent.mt"]] <- PercentageFeatureSet(NMU_O_P, pattern = "^mt-")


# --- Housekeeping gene score ---
# Housekeeping genes are constitutively expressed in all viable cells
# A cell expressing very few housekeeping genes is likely of low quality
# We use a curated mouse housekeeping gene list as an additional QC filter
HK_genes <- read.table(here("data/raw/HK_genes_mouse.txt"))
HK_genes  <- as.vector(HK_genes$V1)

# Intersect the housekeeping list with genes present in our dataset
# (not all housekeeping genes will pass the min.cells filter above)
HK_genes_NMU_O_D <- HK_genes[HK_genes %in% rownames(NMU_O_D)]
HK_genes_NMU_O_P <- HK_genes[HK_genes %in% rownames(NMU_O_P)]

# Count how many housekeeping genes are detected (counts > 0) per cell
# colSums across the HK gene subset gives a per-cell score
NMU_O_D$HK_genes <- colSums(
  GetAssayData(NMU_O_D, assay = "RNA", slot = "counts")[HK_genes_NMU_O_D, ] > 0
)
NMU_O_P$HK_genes <- colSums(
  GetAssayData(NMU_O_P, assay = "RNA", slot = "counts")[HK_genes_NMU_O_P, ] > 0
)


# 3. VISUALIZE QC METRICS

# Violin plots show the distribution of each QC metric across all cells
# We want to see: nFeature_RNA and nCount_RNA normally distributed,
# percent.mt concentrated at low values, HK_genes concentrated at high values
VlnPlot(NMU_O_D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(NMU_O_P, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)

# Scatter plots reveal relationships between QC metrics and help identify:
# - Doublets: cells with abnormally high nCount and nFeature (two cells captured together)
# - Damaged cells: cells with high percent.mt and low HK_genes
plot1 <- FeatureScatter(NMU_O_D, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(NMU_O_D, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot3 <- FeatureScatter(NMU_O_D, feature1 = "percent.mt",  feature2 = "HK_genes")
plot1 + plot2 + plot3

plot1 <- FeatureScatter(NMU_O_P, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(NMU_O_P, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot3 <- FeatureScatter(NMU_O_P, feature1 = "percent.mt",  feature2 = "HK_genes")
plot1 + plot2 + plot3


# 4. FILTER LOW-QUALITY CELLS

# Filtering thresholds applied (based on QC visualization above):
# - nFeature_RNA > 200:  remove empty droplets (too few genes detected)
# - nFeature_RNA < 2500: remove potential doublets (abnormally high gene count)
# - percent.mt < 20:     remove damaged/dying cells
# - HK_genes > 40:       remove cells with insufficient housekeeping gene expression

# Print cell count before filtering to quantify how many cells are removed
NMU_O_D
NMU_O_D <- subset(
  NMU_O_D,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20 & HK_genes > 40
)
NMU_O_D # Print cell count after filtering


NMU_O_P
NMU_O_P <- subset(
  NMU_O_P,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20 & HK_genes > 40
)
NMU_O_P # Print cell count after filtering
