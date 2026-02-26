#✅ 本教程演示了如何使用 AddModuleScore 快速进行代谢通路打分
#✅ 提供了 ViolinPlot、RidgePlot、FeaturePlot、DotPlot、Boxplot、Heatmap 多种可视化方式
#✅ 可直接迁移到 免疫通路、衰老 signature、上皮间质转化 EMT 打分 等其他场景
#AddModuleScore
library(BiocParallel)
library(ggplot2)
library(gtable)
register(MulticoreParam(workers = 8, progressbar = TRUE))
library(Seurat)
#library(SeuratData)
library(tidyverse)
#devtools::install_github("zhanghao-njmu/SCP")
#renv::activate(project = "~/SCP_env")
#library(SCP)
mycolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#FFFF00",
           "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E",
           "#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22",
           "#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
           "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
#data("pbmc3k.final")
setwd("/home/rstudio/project/Metastasis/analysis/Fib_analysis")
load('Pcancer_Fibroblasts_rename_RunUMAP.RData')
scRNA_harmony <- data_harmony
# scRNA_harmony<- UpdateSeuratObject(scRNA_harmony)
scRNA_harmony$celltype <- scRNA_harmony$Subcelltype

# dim(pbmc)
#DimPlot(scRNA_harmony,group.by = 'celltype',cols = mycolour)
# SenMayo_gene <- read.csv("/Users/zhao/Desktop/wechat/衰老基因集/SenMayo_human_senescence_gene.csv",header = T)
# SenMayo_gene <- SenMayo_gene$Gene.human.

# 假设你的 fibroblast 对象叫 fib  #Cancer cell 大综述：Classifying cancer-associated fibroblasts-The good, the bad, and the target
# mCAFs.genes <- c("MMP11", "POSTN", "CTHRC1", "LRRC15", "FAP")
# tCAFs.genes <- c("CD73", "NT5E", "CD10", "MME", "FAP")
# vCAFs.genes <- c("CD146", "MCAM", "SMA", "ACTA2", "NOTCH3", "COL18A1")
# dCAFs.genes <- c("MK167") 
# iCAFs.genes <- c("IL6", "IL1", "CD34", "PLA2G2A", "DPP4", "CFD", "C3", "PI16")
# ifnCAFs.genes <- c("CXCL9", "CXCL10", "CXCL11", "IDO", "ID01")
# rCAFs.genes <- c("CCL19", "CCL2")
# apCAFs.genes <- c("HLA-DR","HLA-DQ","CD74")

mCAFs.genes  <- c("MMP11", "POSTN", "CTHRC1", "LRRC15", "FAP")

tCAFs.genes  <- c("NT5E", "MME", "FAP")

vCAFs.genes  <- c("MCAM", "ACTA2", "NOTCH3", "COL18A1")

dCAFs.genes  <- c("MKI67")

iCAFs.genes  <- c("IL6", "IL1A", "IL1B", "CD34", "PLA2G2A", "DPP4", "CFD", "C3", "PI16")

ifnCAFs.genes <- c("CXCL9", "CXCL10", "CXCL11", "IDO1")

rCAFs.genes <- c("CCL19", "CCL2")

apCAFs.genes <- c("HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "CD74")


scRNA_harmony <- AddModuleScore(object = scRNA_harmony,
                                features = list(mCAFs.genes, tCAFs.genes,vCAFs.genes,dCAFs.genes,iCAFs.genes,ifnCAFs.genes,rCAFs.genes,apCAFs.genes),
                                name = c("mCAFs", "tCAFs","vCAFs","dCAFs","iCAFs","ifnCAFs","rCAFs","apCAFs"))
# 查看打分结果
head(scRNA_harmony@meta.data)

#可视化多种思路
colnames(scRNA_harmony@meta.data)

# 提取 meta.data
df <- scRNA_harmony@meta.data
avg_scores <- df %>%
  group_by(Subcelltype) %>%
  summarise(mCAFs = mean(mCAFs1),
            tCAFs = mean(tCAFs2),
            vCAFs = mean(vCAFs3),
            dCAFs = mean(dCAFs4),
            iCAFs = mean(iCAFs5),
            ifnCAFs = mean(ifnCAFs6),
            rCAFs = mean(rCAFs7),
            apCAFs = mean(apCAFs8)) %>%
  column_to_rownames("Subcelltype")

library(pheatmap)
pheatmap(as.matrix(avg_scores),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Average Pathway Scores by Celltype")
write.csv(avg_scores,file='Fib_celltype_avg_scores.csv')


#1. Violin Plot – 评分分布
VlnPlot(scRNA_harmony,
        features = c("ECM_Score1", "Immunoreg_Score2","Antigen_Score3"),  #"ECM_Score1"
        group.by = "celltype",
        pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#+labs(title = "ECM_Score by Celltype")

# 2. RidgePlot – 分布密度曲线
RidgePlot(scRNA_harmony,
          features = "Lactate_Score1",
          group.by = "celltype") +
  labs(title = "Lactate Score Distribution by Celltype")

#3. FeaturePlot – UMAP 上展示

FeaturePlot(scRNA_harmony,
            features = "Lactate_Score1",
            reduction = "umap",
            cols = c("lightgrey", "red")) +
  labs(title = "Lactate Metabolism Score (UMAP)")

#4. DotPlot – 多通路对比
DotPlot(scRNA_harmony,
        features = c("Lactate_Score1", "FAO_Score2"),
        group.by = "celltype") +
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

#5. Boxplot – 使用 ggplot2 标准箱线图
# 提取 meta.data
df <- scRNA_harmony@meta.data
ggplot(df, aes(x = celltype, y = Lactate_Score1, fill = celltype)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Lactate Metabolism Score by Celltype",
       y = "Lactate_Score1") +
  guides(fill = "none")


#######################Mac
library(BiocParallel)
library(ggplot2)
library(gtable)
register(MulticoreParam(workers = 8, progressbar = TRUE))
library(Seurat)
library(SeuratData)
library(tidyverse)
load('/home/rstudio/project/Metastasis/data/DISCO/0_rawdata/1_anno/norm/finished/decontx/anno/anno02/TOSICA/TOSICA-main/TOSICA_anno/Fib/pancancer_res0.2_anno01_counts_Fib_TOSICAanno_Pro0.3_RunUMAP.RData')


scRNA_harmony <- seu_harmony
rm(seu_harmony)
# scRNA_harmony<- UpdateSeuratObject(scRNA_harmony)
scRNA_harmony$celltype <- scRNA_harmony$anno_Sub
M1 <- c("IL23A","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A","IL1B","IL6",
               "CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7")
M2 <- c("IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24","LYVE1",
              "VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC","CTSD",
              "TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B",
              "FASLG","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4",
              "MRC1","CD163")
ECM <- c("ACTA2","TAGLN","BGN","COL8A1","COL15A1","IGFBP7","TPM1","TPM2","COL10A1",
         "POSTN","MYL9","COL13A1","COL14A1","MYH11","MYLK","ACTG2","FN1","LUM","DCN",
         "VCAN","COL5A1","COL5A2","COL63A","PDGFA","MMP2","COL1A2","RGS5")

Angiogenesis <- c("CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2","FGF18","FGFR1",
                  "FYN","HEY1","ITGAV","JAG1","JAG2","MMP9","NOTCH1","PDGFA","PTK2",
                  "SPP1","STC1","TNFAIP6","TYMP","VAV2","VCAN","VEGFA")

Phagocytosis <-c("MRC1","CD163","MERTK","C1QB") 

Anti_inflammatory <- c("IL1RN","IL10","IL4","IL11","IL13","TGFB1","TNFRSF1A","TNFRSF1B","IL1R2","IL18BP")

Pro_inflammatory <- c("IL1B","TNF","CCL2","CCL3","CCL5","CCL7","CCL8","CCL13","CCL17","CCL22")

IFN <-c("CXCL10","PDL1","ISG15")   #IFN-TAMs（干扰素诱导TAMs）

Reg <-c("ARG1","MRC1","CX3CR1")   #Reg-TAMs（免疫调节型TAMs）

LA <- c("APOC1","APOE","ACP5",'FABP5')     #LA-TAMs（脂质相关TAMs）

RTM <- c("LYVE1","APOE","HES1",'FOLR2')               #RTM-TAMs（组织驻留型TAMs）

Prolif <- c("MKI67")  #Prolif-TAMs（增殖型TAMs）


scRNA_harmony <- AddModuleScore(object = data_harmony,
                                features = list(M1, M2,ECM,Angiogenesis,Phagocytosis,Anti_inflammatory,Pro_inflammatory),
                                name = c("M1_Score", "M2_Score","ECM_Score","Angiogenesis_Score", "Phagocytosis_Score","Anti_inflammatory_Score","Pro_inflammatory_Score"))

# 查看打分结果
head(scRNA_harmony@meta.data)

#可视化多种思路
colnames(scRNA_harmony@meta.data)

###Heatmap – 各 celltype 平均 pathway score

# 提取 meta.data
df <- scRNA_harmony@meta.data
avg_scores <- df %>%
  group_by(celltype) %>%
  summarise(M1_Score = mean(M1_Score1),
            M2_Score = mean(M2_Score2),
            ECM_Score = mean(ECM_Score3),
            Angiogenesis_Score = mean(Angiogenesis_Score4),
            Phagocytosis_Score = mean(Phagocytosis_Score5),
            Anti_inflammatory_Score = mean(Anti_inflammatory_Score6),
            Pro_inflammatory_Score = mean(Pro_inflammatory_Score7)) %>%
  column_to_rownames("celltype")
pheatmap(as.matrix(avg_scores),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Average Pathway Scores by Celltype")
write.csv(avg_scores,file='Mac_celltype_avg_scores.csv')




############others
#1. Violin Plot – 评分分布
VlnPlot(scRNA_harmony,
        features = "Lactate_Score1",
        group.by = "celltype",
        pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Lactate Metabolism Activity by Celltype")
# 2. RidgePlot – 分布密度曲线
RidgePlot(scRNA_harmony,
          features = "Lactate_Score1",
          group.by = "celltype") +
  labs(title = "Lactate Score Distribution by Celltype")
#3. FeaturePlot – UMAP 上展示

FeaturePlot(scRNA_harmony,
            features = "Lactate_Score1",
            reduction = "umap",
            cols = c("lightgrey", "red")) +
  labs(title = "Lactate Metabolism Score (UMAP)")
# 4. DotPlot – 多通路对比
DotPlot(scRNA_harmony,
        features = c("Lactate_Score1", "FAO_Score2"),
        group.by = "celltype") +
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
#5. Boxplot – 使用 ggplot2 标准箱线图
# 提取 meta.data
df <- scRNA_harmony@meta.data
ggplot(df, aes(x = celltype, y = Lactate_Score1, fill = celltype)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Lactate Metabolism Score by Celltype",
       y = "Lactate_Score1") +
  guides(fill = "none")


####################################CIN70 
# Some required R packages
library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(cowplot)
library(ggplot2)

# Parallel Computing
library(future)
availableWorkers()
availableCores()
options(future.globals.maxSize= 30 * 1024^3)
plan(strategy = "multicore", workers = 30)

# Step1 : Define gene set
CIN70 <- c(
  "TPX2", "PRC1", "FOXM1", "CDK1", "RAB5IF", "TGIF2", "MCM2", "H2AFZ", "TOP2A",
  "PCNA", "UBE2C", "MELK", "TRIP13", "NCAPD2", "MCM7", "RNASEH2A", "RAD51AP1", "KIF20A",
  "CDC45", "MAD2L1", "ESPL1", "CCNB2", "FEN1", "TTK", "CCT5", "RFC4", "ATAD2", "CKAP5",
  "NUP205", "CDC20", "CKS2", "RRM2", "ELAVL1", "CCNB1", "RRM1", "AURKB", "EZH2", "CTPS1",
  "DKC1", "OIP5", "CDCA8", "PTTG1", "CEP55", "H2AFX", "CMAS", "NCAPH", "MCM10", "LSM4",
  "NCAPG2", "ASF1B", "ZWINT", "PBK", "ZWILCH", "CDCA3", "ECT2", "CDC6", "UNG", "MTCH2",
  "RAD21", "ACTL6A", "GPI", "PDCD2L", "SRSF2", "HDGF", "NXT1", "NEK2", "DHCR7", "AURKA",
  "NDUFAB1", "NEMP1", "KIF4A"
)

# Step2 : AddModuleScore
geneset <- list("CIN70"=CIN70)

scObject <- AddModuleScore(object = scObject, features = geneset, name = names(geneset))

# Step3 : Plots
v1 <- VlnPlot(object = scObject, features = "CIN70", group.by="Disease", pt.size=0)+
  geom_boxplot(outlier.size=0, color="black", width=0.3)

v2 <- VlnPlot(object = scObject, features = "CIN70", group.by="TumorGroup_Disease", pt.size=0)+
  geom_boxplot(outlier.size=0, color="black", width=0.3)


