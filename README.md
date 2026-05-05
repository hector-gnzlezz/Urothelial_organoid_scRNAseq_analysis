# Urothelial Organoid scRNA-seq Analysis

Single-cell RNA-seq analysis of mouse urothelial organoids comparing 
differentiated (NMU-O-D) and proliferating (NMU-O-P) conditions to 
characterize cell populations and study the urothelial differentiation program.

## Biological context

Urothelial organoids recapitulate the architecture of bladder epithelium 
in vitro. This analysis identifies and characterizes basal, intermediate 
and luminal cell populations, and quantifies how their proportions shift 
between differentiated and proliferating culture conditions.

Single-cell findings are projected onto bulk RNA-seq data via:
- **ssGSEA** — scoring bulk samples for cell population activity
- **CIBERSORTx** — estimating cell type fractions in bulk samples

## Data

Publicly available at GEO: [GSE131909](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131909)  
Publication: [PMID 31562298](https://pubmed.ncbi.nlm.nih.gov/31562298/)  
See `data/README.md` for download instructions.

## Pipeline

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_QC_filtering.R` | Load 10X data, compute QC metrics, filter low-quality cells |
| 2 | `02_normalization.R` | Normalize, find variable features, PCA, UMAP, cell cycle scoring, doublet detection |
| 3 | `03_clustering.R` | Resolution sweep, marker identification, cell type annotation |
| 4 | `04_sample_integration.R` | Harmony batch correction, integrated clustering, gene signature export |

## Key results

- Identified 4 cell populations: Basal, Intermediate, Luminal, Basal G2M
- Differentiated condition enriched in Luminal and Intermediate cells
- Proliferating condition shows expansion of Basal G2M population
- Cell population gene signatures exported for ssGSEA and CIBERSORTx

## How to run

```r
# 1. Restore R environment
install.packages("renv")
renv::restore()

# 2. Download raw data (see data/README.md)

# 3. Run pipeline in order
source("scripts/01_QC_filtering.R")
source("scripts/02_normalization.R")
source("scripts/03_clustering.R")
source("scripts/04_sample_integration.R")
```

## Repository structure
```
├── data/
│   ├── raw/                  # 10X matrices (tracked; also available at GEO GSE131909)
│   │   └── HK_genes_mouse.txt
│   └── processed/            # Pipeline outputs (not tracked)
├── results/
│   └── figures/              # Generated plots
├── scripts/
│   ├── 01_QC_filtering.R
│   ├── 02_normalization.R
│   ├── 03_clustering.R
│   └── 04_sample_integration.R
└── renv.lock                 # Reproducible R environment
```
## Dependencies

R 4.x — managed with `renv`. Key packages: Seurat, harmony, DoubletFinder, 
nichenetr, clustree.
