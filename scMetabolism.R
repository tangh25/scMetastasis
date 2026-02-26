
# devtools::install_github("YosefLab/VISION")
# devtools::install_github("wu-yc/scMetabolism")
# BiocManager::install('AUCell')
# BiocManager::install('GSEABase')
# BiocManager::install('GSVA')
library(GSVA)

library(scMetabolism)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(Seurat)
library(tidyverse)
library(BiocParallel)
register(MulticoreParam(workers = 8, progressbar = TRUE)) 
load('RunUMAP.RData')    #sce.all
Idents(sce.all) <- "Sub_celltype"  ##celltype
set.seed(20250103)  
countexp.Seurat <- sc.metabolism.Seurat(obj = sce.all, method = "AUCell", imputation = F, ncores = 2, metabolism.type = "KEGG")
input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[1:20]
pdf("./scMetabolism_Fib.pdf",height = 10,width =10)
scMetabolism::DotPlot.metabolism(obj = countexp.Seurat,                    
                                 pathway = input.pathway,                    
                                 phenotype = "Subcelltype",                    
                                 norm = "y")+
  theme(axis.text = element_text(color = "black"))
dev.off()
saveRDS(countexp.Seurat,file='RunUMAP_scMetabolism.RData')
