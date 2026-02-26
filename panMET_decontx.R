##decontx
# Clear workspace
rm(list = ls())
# Load required libraries
library(Seurat)     # Single-cell analysis framework
library(decontX)    # Ambient RNA contamination removal
library(dplyr)      # Data manipulation

# Set working directory containing input .rds files
setwd("/home/rstudio/project/Metastasis/data/DISCO/0_rawdata/1_anno/GEO25")

# Get all .rds files with "norm.rds" suffix
rds_files <- list.files(pattern = "norm\\.rds$")

# Loop through each RDS file
for (file in rds_files) {
  
  cat("Processing file:", file, "\n")
  
  # Read RDS file (assumed to contain count matrix)
  sce1 <- readRDS(file)
  
  # Ensure gene names (rownames) exist
  if (is.null(rownames(sce1))) {
    stop(paste("File", file, "has missing rownames. Please check data format."))
  }

  # Create Seurat object from count matrix
  sce <- CreateSeuratObject(
    counts = sce1, 
    min.cells = 3,       # Keep genes expressed in ≥3 cells
    min.features = 200   # Keep cells with ≥200 detected genes
  )
  
  # Extract raw counts matrix (Seurat v4 format)
  set.seed(123)
  counts <- sce@assays$RNA@counts
  
  # Run decontX to estimate ambient RNA contamination
  decontX_results <- decontX(counts)
  
  # Add contamination score to metadata
  sce$contamination <- decontX_results$contamination
  
  # Filter out cells with high contamination (≥ 20%)
  sce_filt <- sce[, sce$contamination < 0.2]
  
  # Generate new output filename
  new_file_name <- sub("\\.rds$", "_rnafilt.rds", file)
  
  # Save filtered Seurat object
  saveRDS(sce_filt, file = new_file_name)
  
  cat("Finished:", new_file_name, "\n")
}

# Final message after all files processed
cat("All files processed successfully.\n")
