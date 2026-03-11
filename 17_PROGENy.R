########################## PROGENy Pathway Activity Analysis #################################
# Load packages
library(Seurat)
library(tidyverse)
library(data.table)
library(progeny)
library(pheatmap)
library(ggpubr)
library(ggsci)
library(viridis)

# Set working directory and load Seurat object
setwd("./analysis")
load('data_RunUMAP.RData')

# Preprocess
sc_data$celltype <- sc_data$Subcelltype
sce <- sc_data
rm(sc_data)
Idents(sce) <- "celltype"

# Create cluster annotation dataframe
CellsClusters <- data.frame(Cell = names(Idents(sce)), 
                            CellType = as.character(Idents(sce)),
                            stringsAsFactors = FALSE)

# UMAP visualization
DimPlot(sce, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Run PROGENy
sce <- progeny(sce, scale = FALSE, organism = "Human", top = 500, perm = 1, 
               return_assay = TRUE)  # Results saved in a separate assay
sce <- Seurat::ScaleData(sce, assay = "progeny") 

# Extract pathway activity scores
progeny_scores_df <- as.data.frame(t(GetAssayData(sce, slot = "scale.data", assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

# Join with cell type information
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

# Summarize scores by pathway and cell type
summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

# Reshape for heatmap plotting
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# Save summarized scores
write.csv(summarized_progeny_scores_df, file = 'summarized_progeny_scores_df.csv')

########################## Heatmap #################################
paletteLength <- 100
myColor <- colorRampPalette(c("#5e8ac2", "white","#c04343"))(paletteLength)
progenyBreaks <- c(
  seq(min(summarized_progeny_scores_df), 0, length.out = ceiling(paletteLength/2) + 1),
  seq(max(summarized_progeny_scores_df)/paletteLength, max(summarized_progeny_scores_df),
      length.out = floor(paletteLength/2))
)
pdf("fib_progeny2.pdf", height = 5, width = 8)
pheatmap(t(summarized_progeny_scores_df[,-1]), fontsize = 14, fontsize_row = 10,
         color = myColor, breaks = progenyBreaks,
         main = "PROGENy (500)", angle_col = 45,
         treeheight_col = 0, border_color = NA)
dev.off()
