##http://172.25.51.2:8716/   ###sc.metabolism
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
setwd("/home/rstudio/project/Metastasis/data/DISCO/0_rawdata/1_anno/norm/finished/decontx/anno/anno02/adata_by_celltype/")
setwd('/home/rstudio/project/Metastasis/analysis/Cottrazm/TumorBoundary/bin')
#load("./data_harmony_final.RData")
sce.all = schard::h5ad2seurat('Fibroblasts.h5ad')
saveRDS(sce.all,file = "Fibroblasts.rds")
sce.all=readRDS("pancancer_res0.2_anno01_counts_Fib_TOSICAanno_Pro0.3_RunUMAP.RDataPcancer_anno03_sub_celltype_anno_clinic1_MET_Fib.rds")
load('Macrophages1_RunUMAP.RData')
#清除无关细胞类型
sce.all <- sce.all[, !sce.all$anno_sub %in% c("Mesenchymal stem-like cells", "SMC", "Pericyte")]
sce.all$anno_sub[sce.all$anno_sub == "myCAFs"] <- "SOD2+ myCAFs"
Idents(sce.all) <- "anno_Sub"  ##celltype
data_sub <- sce.all
rm(sce.all)
# 1. 标准预处理流程
set.seed(20250413)
data_sub <- NormalizeData(data_sub,normalization.method = "LogNormalize", scale.factor = 10000)
data_sub <- FindVariableFeatures(data_sub,selection.method = "vst", nfeatures = 3000)
data_sub <- ScaleData(data_sub,features = rownames(data_sub))
data_sub <- Seurat::RunPCA(data_sub, features = VariableFeatures(object = data_sub))

set.seed(20250413)
library(harmony)
dat_harmony <- data_sub %>% RunHarmony("sample", plot_convergence = T,max_iter=20,reduction.save ="harmony")
ElbowPlot(dat_harmony,ndims=50, reduction="harmony")

data_harmony <- FindNeighbors(dat_harmony, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.5)
data_harmony <- RunUMAP(data_harmony, reduction = "harmony", dims = 1:30)
data_harmony <- RunTSNE(data_harmony, reduction = "harmony", dims = 1:30)
save(data_harmony,file ="./Fib_RunUMAP.RData")
UMAP_celltype <- DimPlot(data_harmony, reduction ="umap",
                         group.by="anno_Sub",label = F);UMAP_celltype
Idents(data_harmony) <- data_harmony$anno_Sub
scCustomize::DimPlot_scCustom(data_harmony, figure_plot = TRUE)
clusterCols <- c(
  '#E5D2DD','#53A85F','#F1BB72','#F3B1A0','#D6E7A3',
  '#57C3F3','#E95C59','#E59CC4','#AB3282','#BD956A',
  '#9FA3A8','#E0D4CA','#C5DEBA','#F7F398','#C1E6F3',
  '#6778AE','#91D0BE','#B53E2B','#712820','#DCC1DD',
  '#CCE0F5','#CCC9E6','#625D9E','#68A180','#3A6963',
  '#968175','#F0E68C','#FFFFE0','#EE82EE','#FF6347',
  '#1E90FF','#20B2AA','#FF4500','#DA70D6','#B0C4DE'
)

p <- DimPlot(data_harmony, reduction = "umap",label = T,
             raster = F,cols = clusterCols)+
  theme_bw()+
  theme(plot.title = element_text(size=15,hjust=0.5),
        legend.position = "right",
        legend.title = element_text(colour="black", size=10,face="plain"),
        legend.text = element_text(colour="black", size=9.5,face ="plain"),
        legend.background = element_rect(fill="white"),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))+
  guides(color = guide_legend(override.aes = list(size = 3),ncol = 2,byrow = F,reverse = F))
ggsave("Fib_umap.pdf", plot = p, width = 6, height = 4)


###sc.metabolism  PMID: 39367943 Figure4I
rm(dat_harmony)
rm(data_sub)
set.seed(20250103)
countexp.Seurat <- sc.metabolism.Seurat(obj = data_harmony, method = "AUCell", imputation = F, ncores = 2, metabolism.type = "KEGG")
input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[1:20]
pdf("./scMetabolism_Fib.pdf",height = 10,width =10)
scMetabolism::DotPlot.metabolism(obj = countexp.Seurat,                    
                                 pathway = input.pathway,                    
                                 phenotype = "Subcelltype",                    
                                 norm = "y")+
  theme(axis.text = element_text(color = "black"))
dev.off()
saveRDS(countexp.Seurat,file='Fib_RunUMAP_scMetabolism.RData')
################################################################################

