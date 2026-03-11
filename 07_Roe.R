# Load required libraries
library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)
library(Seurat)
library(tidyverse)
library(readr)
#library(qs)
library(BiocParallel)
register(MulticoreParam(workers = 10, progressbar = TRUE)) 
setwd("")

setwd("./analysis")
############################
# UMAP visualization
############################

DimPlot(sc_data, reduction = "umap")
DimPlot(sc_data, reduction = "umap", group.by = "Subcelltype")  # visualize by subcell type

# Inspect object structure
str(sc_data@reductions$umap)
str(sc_data@meta.data)
table(Idents(sc_data))

DimPlot(sc_data)

############################
# Export metadata
############################

data  <- sc_data@meta.data
colnames(data)

write.csv(Roe,file = 'roe_metastasis_tissue.csv')

setwd("./analysis")

data3 <- read.csv('Pcance_rename_metadata.csv')


############################
# Prepare metadata for Ro/e analysis
# Required columns:
# clone.id (cell ID), patient (sample), majorCluster (cell subtype), loc (tissue/location)
############################

#### Metastasis vs Primary comparison

data <- data[,c(4,5,33,19)]
data <- data3[,c(5,6,30,19)]

colnames(data) <- c("cell_id","sample","anno_major","group")

data$Cell_Name <- c("Major cell")

data$patient <- as.character(data$sample)
data$clone.id <- as.character(data$cell_id)
data$majorCluster <- as.character(data$anno_major)
data$loc <- as.character(data$group)


############################
# Tissue-specific analysis
############################

data  <- sc_data@meta.data
colnames(data)

data <- data[,c(4,5,32,19)]

colnames(data) <- c("cell_id","sample","Subcelltype","tissue")

data$Cell_Name <- c("Mac cell")

data$patient <- as.character(data$sample)
data$clone.id <- as.character(data$cell_id)
data$majorCluster <- as.character(data$Subcelltype)
data$loc <- as.character(data$tissue)


############################
# Subset metastasis samples
############################

data  <- sc_data@meta.data

data1 <- subset(data3,group=='Metastasis')

colnames(data1)

data1 <- data1[,c(5,6,34,32)]

colnames(data1) <- c("cell_id","sample","Subcelltype","tissue")

data1$Cell_Name <- c("Mac cell")

data1$patient <- as.character(data1$sample)
data1$clone.id <- as.character(data1$cell_id)
data1$majorCluster <- as.character(data1$Subcelltype)
data1$loc <- as.character(data1$tissue)


############################
# Ro/e analysis
############################

##### Metastasis vs Primary

Roe <- calTissueDist(
                     data1,
                     byPatient = F,
                     colname.cluster = "majorCluster", # cell subtype
                     colname.patient = "patient",      # sample
                     colname.tissue = "loc",           # tissue/location
                     method = "chisq",                 # chisq, fisher, or freq
                     min.rowSum = 0
                     )

Roe

write.csv(Roe,file = 'roe_metastasis_groupOM.csv')

table(data1$majorCluster)


############################
# Heatmap visualization (numeric version)
############################

col_fun <- colorRamp2(
                      c(min(Roe, na.rm = TRUE), 1, max(Roe, na.rm = TRUE)), 
                      c("blue", "white", "red")
                      )

Heatmap(as.matrix(Roe),
        show_heatmap_legend = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        row_names_side = 'right', 
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "Ro/e Index",
          at = seq(0.5, 1.5, by = 0.5),
          labels = seq(0.5, 1.5, by = 0.5)
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", Roe[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
        }
)


############################
# Heatmap visualization (symbol version)
############################

col_fun  <-  colorRamp2(
                        c(min(Roe,na.rm  =  TRUE),1,max(Roe,na.rm  =TRUE)),  
                        c("#f6f8e6", "#f9a33e", "red")
                        )

Heatmap(as.matrix(Roe),
        show_heatmap_legend = TRUE, 
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = 'right',
        column_names_side = "top",
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12),
        heatmap_legend_param = list(
          title = "Ro/e",
          at = c(0, max(Roe)),
          labels = c("0", "Max.")
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {

          value <- Roe[i, j]

          # Symbol annotation based on Ro/e value
          symbol <- if(value == 0) {
            "−"
          } else if(value > 0 & value < 0.2) {
            "+/−"
          } else if(value >= 0.2 & value <= 0.8) {
            "+"
          } else if(value > 0.8 & value <= 1) {
            "++"
          } else if(value > 1) {
            "+++"
          }

          grid.text(symbol, x, y, gp = gpar(fontsize = 10, col = "black"))
        }
)


############################
# Save results
############################

Roe=as.data.frame(Roe)

write_csv(Roe,file='Roe_met_tissue.csv')

saveRDS(data,file='Roe_data.rds')
