# ==============================
# irGSEA installation section
# ==============================

# Install AUCell from CRAN
install.packages('AUCell')

# List of CRAN packages required for irGSEA workflow
cran.packages <- c("aplot", "BiocManager", "data.table", "devtools", 
                   "doParallel", "doRNG", "dplyr", "ggfun", "gghalves", 
                   "ggplot2", "ggplotify", "ggridges", "ggsci", "irlba",
                   "magrittr", "Matrix", "msigdbr", "pagoda2", "pointr", 
                   "purrr", "RcppML", "readr", "reshape2", "reticulate", 
                   "rlang", "RMTstat", "RobustRankAggreg", "roxygen2", 
                   "Seurat", "SeuratObject", "stringr", "tibble", "tidyr", 
                   "tidyselect", "tidytree", "VAM")

# Install missing CRAN packages
for (i in cran.packages) {
  if (!requireNamespace(i, quietly = TRUE)) {
    install.packages(i, ask = FALSE, update = FALSE)
  }
}

# List of Bioconductor packages
bioconductor.packages <- c("AUCell", "BiocParallel", "ComplexHeatmap", 
                           "decoupleR", "fgsea", "ggtree", "GSEABase", 
                           "GSVA", "Nebulosa", "scde", "singscore",
                           "SummarizedExperiment", "UCell",
                           "viper","sparseMatrixStats")

# Install missing Bioconductor packages
for (i in bioconductor.packages) {
  if (!requireNamespace(i, quietly = TRUE)) {
    install.packages(i, ask = FALSE, update = FALSE)
  }
}

# Install irGSEA from GitHub if not installed
if (!requireNamespace("irGSEA", quietly = TRUE)) { 
  devtools::install_github("chuiqin/irGSEA", force = TRUE)
}

########################################
# irGSEA analysis workflow
########################################

rm(list = ls())

library(irGSEA)
library(Seurat)
library(RcppML)
library(scCustomize)
library(doMC)

# Set working directory
setwd("/home/rstudio/project/Metastasis/analysis/Fib_analysis")

# Load preprocessed fibroblast Seurat object
load('RunUMAP.RData')  #scRNA

# Set celltype annotation
scRNA$celltype <- scRNA$Sub_celltype
Idents(scRNA) <- "celltype"

# Assign to working object
sc_dataset <- scRNA
rm(scRNA)

# Visualize UMAP colored by cell type
UMAP_celltype <- DimPlot(sc_dataset, reduction ="umap",
                         group.by="celltype", label = FALSE)
UMAP_celltype

# Ensure identities are set correctly
Idents(sc_dataset) <- sc_dataset$celltype

########################################################
# irGSEA enrichment score calculation
########################################################

# Update to Seurat v5 object format
sc_dataset <- SeuratObject::UpdateSeuratObject(sc_dataset)

# Normalize expression data
sc_dataset2 <- NormalizeData(sc_dataset)

rm(sc_dataset)

# Compute enrichment scores using multiple scoring methods
sc_dataset2 <- irGSEA.score(
  object = sc_dataset2,
  assay = "RNA",
  slot = "data",
  seeds = 123,
  min.cells = 3,
  min.feature = 0,
  custom = FALSE,
  geneset = NULL,
  msigdb = TRUE,
  species = "Homo sapiens",
  category = "H",      # MSigDB Hallmark gene sets
  subcategory = NULL,
  geneid = "symbol",
  method = c("AUCell","UCell","singscore",
             "ssgsea", "JASMINE", "viper"),
  aucell.MaxRank = NULL,
  ucell.MaxRank = NULL,
  kcdf = 'Gaussian'
)

########################################################
# Integrate differential gene sets across scoring methods
########################################################

# Increase memory limit if necessary
options(future.globals.maxSize = 100000 * 1024^5)

# Perform differential enrichment analysis across cell types
result.dge <- irGSEA.integrate(
  object = sc_dataset2,
  group.by = "celltype",
  method = c("AUCell","UCell","singscore",
             "ssgsea", "JASMINE", "viper")
)

# Extract gene sets identified by RRA (Robust Rank Aggregation)
geneset.show <- result.dge$RRA %>% 
  dplyr::filter(pvalue <= 0.05) %>% 
  dplyr::pull(Name) %>% 
  unique()

########################################################
# Visualization
########################################################

# Heatmap of top 10 gene sets based on RRA results
heatmap.plot <- irGSEA.heatmap(
  object = result.dge,
  method = "RRA",
  top = 10,
  show.geneset = NULL
)
heatmap.plot

# Save integrated differential enrichment results
saveRDS(result.dge, file='result.dge_irGSEA_RRA.rds')

# Heatmap using ssgsea scoring method
heatmap.plot <- irGSEA.heatmap(
  object = result.dge,
  method = "ssgsea",
  top = 10,
  show.geneset = NULL
)
