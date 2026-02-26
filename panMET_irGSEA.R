#http://172.25.51.2:8716/
#irGSEA
# install packages from CRAN
install.packages('AUCell')
cran.packages <- c("aplot", "BiocManager", "data.table", "devtools", 
                   "doParallel", "doRNG", "dplyr", "ggfun", "gghalves", 
                   "ggplot2", "ggplotify", "ggridges", "ggsci", "irlba",
                   "magrittr", "Matrix", "msigdbr", "pagoda2", "pointr", 
                   "purrr", "RcppML", "readr", "reshape2", "reticulate", 
                   "rlang", "RMTstat", "RobustRankAggreg", "roxygen2", 
                   "Seurat", "SeuratObject", "stringr", "tibble", "tidyr", 
                   "tidyselect", "tidytree", "VAM")

for (i in cran.packages) {
  if (!requireNamespace(i, quietly = TRUE)) {
    install.packages(i, ask = F, update = F)
  }
}

# install packages from Bioconductor
bioconductor.packages <- c("AUCell", "BiocParallel", "ComplexHeatmap", 
                           "decoupleR", "fgsea", "ggtree", "GSEABase", 
                           "GSVA", "Nebulosa", "scde", "singscore",
                           "SummarizedExperiment", "UCell",
                           "viper","sparseMatrixStats")

for (i in bioconductor.packages) {
  if (!requireNamespace(i, quietly = TRUE)) {
    install.packages(i, ask = F, update = F)
  }
}

# install packages from Github
if (!requireNamespace("irGSEA", quietly = TRUE)) { 
  devtools::install_github("chuiqin/irGSEA", force =T)
}

########################################

rm(list = ls())
library(irGSEA)
library(Seurat)
#library(SeuratData)
library(RcppML)
library(scCustomize)
#BiocManager::install('doMC')
library('doMC')
setwd("/home/rstudio/project/Metastasis/analysis/Fib_analysis")
load('Pcancer_Fibroblasts_rename_RunUMAP.RData')
#sc_dataset <- readRDS("./pancancer_res0.2_anno01_counts_Fib_TOSICAanno_Pro0.3.rds")
data_harmony$celltype <- data_harmony$Subcelltype
Idents(data_harmony) <- "celltype"
sc_dataset <-data_harmony
rm(data_harmony)
# sc_dataset <- scRNA
# Check
UMAP_celltype <- DimPlot(sc_dataset, reduction ="umap",
                         group.by="celltype",label = F);UMAP_celltype
Idents(sc_dataset) <- sc_dataset$celltype
#scCustomize::DimPlot_scCustom(sc_dataset, figure_plot = TRUE)


#irGSEA计算富集分数
#### Seurat V5对象 ####
sc_dataset <- SeuratObject::UpdateSeuratObject(sc_dataset)
# sc_dataset2 <- CreateSeuratObject(counts = CreateAssay5Object(GetAssayData(sc_dataset,
#                                                                            assay = "RNA", 
#                                                                            slot="counts")),
#                                   meta.data = sc_dataset[[]])
#sc_dataset2 <- sc_dataset
sc_dataset2 <- NormalizeData(sc_dataset)
rm(sc_dataset)
rm(data_harmony)
sc_dataset2 <- irGSEA.score(object = sc_dataset2, assay = "RNA",
                            slot = "data", seeds = 123, 
                            #ncores = 1,
                            min.cells = 3, min.feature = 0,
                            custom = F, geneset = NULL, msigdb = T,
                            species = "Homo sapiens", 
                            category = "H",  
                            subcategory = NULL, 
                            geneid = "symbol",
                            method = c("AUCell","UCell","singscore",
                                       "ssgsea", "JASMINE", "viper"),
                            aucell.MaxRank = NULL, 
                            ucell.MaxRank = NULL,
                            kcdf = 'Gaussian')


# 整合差异基因集
# 如果报错，考虑加句代码
options(future.globals.maxSize = 100000 * 1024^5)
result.dge <- irGSEA.integrate(object = sc_dataset2,
                               group.by = "celltype",
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea", "JASMINE", "viper"))

# 查看RRA识别的在多种打分方法中都普遍认可的差异基因集
geneset.show <- result.dge$RRA %>% 
  dplyr::filter(pvalue <= 0.05) %>% 
  dplyr::pull(Name) %>% unique(.)

#4、可视化
#热图
heatmap.plot <- irGSEA.heatmap(object = result.dge,
                               method = "RRA", 
                               top = 10,
                               show.geneset = NULL)
heatmap.plot
saveRDS(result.dge,file='fib_result.dge_irGSEA_RRA.rds')
#
heatmap.plot <- irGSEA.heatmap(object = result.dge,
                               method = "ssgsea", 
                               top = 10,
                               show.geneset = NULL)

heatmap.plot

#####################
heatmap.plot <- irGSEA.heatmap(object = result.dge,
                              method = "RRA", #从'RRA"换成“ssgsea”
                              top = 10,
                              show.geneset = geneset.show)
heatmap.plot

#气泡图
bubble.plot <- irGSEA.bubble(object = result.dge,
                             method = "RRA",
                             top = 10,
                             show.geneset = geneset.show)
bubble.plot

#Upset图
upset.plot <- irGSEA.upset(object = result.dge,
                           method = "RRA")
upset.plot
#堆叠条形图
barplot.plot <- irGSEA.barplot(object = result.dge,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea", "JASMINE", "viper"))
barplot.plot

Idents(sc_dataset2) <- sc_dataset2$celltype
halfvlnplot <- irGSEA.halfvlnplot(object = sc_dataset2,
                                  method = "AUCell",
                                  show.geneset = "HALLMARK-NOTCH-SIGNALING")
halfvlnplot

vlnplot <- irGSEA.vlnplot(object = sc_dataset2,
                          method = c("AUCell", "UCell", 
                                     "singscore", "ssgsea", 
                                     "JASMINE", "viper"),
                          show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE")
vlnplot
#山峦图
ridgeplot <- irGSEA.ridgeplot(object = sc_dataset2,
                              method = "AUCell",
                              show.geneset = "HALLMARK-NOTCH-SIGNALING")
ridgeplot

#密度热图
densityheatmap <- irGSEA.densityheatmap(object = sc_dataset2,
                                        method = "AUCell",
                                        show.geneset = "HALLMARK-NOTCH-SIGNALING")
densityheatmap
