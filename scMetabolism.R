# Install required packages if not already installed
# devtools::install_github("YosefLab/VISION")
# devtools::install_github("wu-yc/scMetabolism")
# BiocManager::install('AUCell')
# BiocManager::install('GSEABase')
# BiocManager::install('GSVA')

library(GSVA)            # Gene set variation analysis methods
library(scMetabolism)    # Single-cell metabolism analysis toolkit
library(ggplot2)         # Visualization
library(ggpubr)          # Publication-ready plots
library(ComplexHeatmap)  # Advanced heatmaps
library(RColorBrewer)    # Color palettes
library(tidyverse)       # Data manipulation
library(Seurat)          # Single-cell analysis framework
library(BiocParallel)    # Parallel computing

# Register parallel backend (8 CPU cores)
register(MulticoreParam(workers = 8, progressbar = TRUE)) 

# Load preprocessed Seurat object with UMAP embedding
load('RunUMAP.RData')    # sce.all

# Set cell identity to Sub_celltype annotation
Idents(sce.all) <- "Sub_celltype"

# Set random seed for reproducibility
set.seed(20250103)  

# Calculate metabolic pathway activity scores using AUCell method
# KEGG metabolism gene sets are used
countexp.Seurat <- sc.metabolism.Seurat(
  obj = sce.all,
  method = "AUCell",
  imputation = FALSE,
  ncores = 2,
  metabolism.type = "KEGG"
)

# Select the first 20 metabolic pathways for visualization
input.pathway <- rownames(
  countexp.Seurat@assays[["METABOLISM"]][["score"]]
)[1:20]

# Open PDF device for saving the plot
pdf("./scMetabolism.pdf", height = 10, width = 10)

# Generate dot plot of metabolic activity across cell subtypes
scMetabolism::DotPlot.metabolism(
  obj = countexp.Seurat,
  pathway = input.pathway,
  phenotype = "Subcelltype",   # Cell grouping variable
  norm = "y"                   # Normalize across pathways
) +
  theme(axis.text = element_text(color = "black"))

# Close PDF device
dev.off()

# Save the metabolism-scored Seurat object for downstream analysis
saveRDS(countexp.Seurat, file = 'RunUMAP_scMetabolism.RData')
