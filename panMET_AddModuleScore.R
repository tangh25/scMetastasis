# AddModuleScore for CAF-related gene signatures

library(BiocParallel)
library(ggplot2)
library(gtable)

# Register parallel computing (8 cores)
register(MulticoreParam(workers = 8, progressbar = TRUE))

library(Seurat)
library(tidyverse)

# Custom color palette (not used later but kept for plotting)
mycolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#FFFF00",
           "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E",
           "#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22",
           "#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
           "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

# Load preprocessed Seurat object with UMAP results
load('RunUMAP.RData')  # scRNA object

# Use Sub_celltype as celltype annotation
scRNA$celltype <- scRNA$Sub_celltype

# CAF subtype gene signatures
mCAFs.genes  <- c("MMP11", "POSTN", "CTHRC1", "LRRC15", "FAP")
tCAFs.genes  <- c("NT5E", "MME", "FAP")
vCAFs.genes  <- c("MCAM", "ACTA2", "NOTCH3", "COL18A1")
dCAFs.genes  <- c("MKI67")
iCAFs.genes  <- c("IL6", "IL1A", "IL1B", "CD34", "PLA2G2A", "DPP4", "CFD", "C3", "PI16")
ifnCAFs.genes <- c("CXCL9", "CXCL10", "CXCL11", "IDO1")
rCAFs.genes <- c("CCL19", "CCL2")
apCAFs.genes <- c("HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "CD74")

# Compute module scores for each CAF subtype
scRNA <- AddModuleScore(
  object = scRNA,
  features = list(mCAFs.genes, tCAFs.genes, vCAFs.genes, dCAFs.genes,
                  iCAFs.genes, ifnCAFs.genes, rCAFs.genes, apCAFs.genes),
  name = c("mCAFs", "tCAFs", "vCAFs", "dCAFs",
           "iCAFs", "ifnCAFs", "rCAFs", "apCAFs")
)

# Inspect metadata with scores
head(scRNA@meta.data)
colnames(scRNA@meta.data)

# Extract metadata
df <- scRNA@meta.data

# Calculate average module scores per fibroblast subtype
avg_scores <- df %>%
  group_by(Sub_celltype) %>%
  summarise(
    mCAFs = mean(mCAFs1),
    tCAFs = mean(tCAFs2),
    vCAFs = mean(vCAFs3),
    dCAFs = mean(dCAFs4),
    iCAFs = mean(iCAFs5),
    ifnCAFs = mean(ifnCAFs6),
    rCAFs = mean(rCAFs7),
    apCAFs = mean(apCAFs8)
  ) %>%
  column_to_rownames("Sub_celltype")

library(pheatmap)

# Plot heatmap of scaled average scores
pheatmap(
  as.matrix(avg_scores),
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Average Pathway Scores by Celltype"
)

# Save results to CSV
write.csv(avg_scores, file = 'Fib_celltype_avg_scores.csv')


####################### Macrophage signatures

library(BiocParallel)
library(ggplot2)
library(gtable)

# Register parallel computing
register(MulticoreParam(workers = 8, progressbar = TRUE))

library(Seurat)
library(SeuratData)
library(tidyverse)

# Load Seurat object
load('RunUMAP.RData')  # scRNA object

# Use Sub_celltype as celltype annotation
scRNA$celltype <- scRNA$Sub_celltype

# Define macrophage functional gene signatures

M1 <- c("IL23A","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A","IL1B","IL6",
        "CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7")

M2 <- c("IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24","LYVE1",
        "VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC","CTSD",
        "TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B",
        "FASLG","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4",
        "MRC1","CD163")

ECM <- c("ACTA2","TAGLN","BGN","COL8A1","COL15A1","IGFBP7","TPM1","TPM2",
         "COL10A1","POSTN","MYL9","COL13A1","COL14A1","MYH11","MYLK",
         "ACTG2","FN1","LUM","DCN","VCAN","COL5A1","COL5A2","COL63A",
         "PDGFA","MMP2","COL1A2","RGS5")

Angiogenesis <- c("CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2",
                  "FGF18","FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2",
                  "MMP9","NOTCH1","PDGFA","PTK2","SPP1","STC1","TNFAIP6",
                  "TYMP","VAV2","VCAN","VEGFA")

Phagocytosis <- c("MRC1","CD163","MERTK","C1QB")

Anti_inflammatory <- c("IL1RN","IL10","IL4","IL11","IL13","TGFB1",
                       "TNFRSF1A","TNFRSF1B","IL1R2","IL18BP")

Pro_inflammatory <- c("IL1B","TNF","CCL2","CCL3","CCL5","CCL7",
                      "CCL8","CCL13","CCL17","CCL22")

# Additional TAM-related gene sets (not used in scoring below)
IFN <- c("CXCL10","PDL1","ISG15")
Reg <- c("ARG1","MRC1","CX3CR1")
LA <- c("APOC1","APOE","ACP5",'FABP5')
RTM <- c("LYVE1","APOE","HES1",'FOLR2')
Prolif <- c("MKI67")

# Compute module scores for macrophage functions
scRNA <- AddModuleScore(
  object = scRNA,
  features = list(M1, M2, ECM, Angiogenesis,
                  Phagocytosis, Anti_inflammatory, Pro_inflammatory),
  name = c("M1_Score", "M2_Score", "ECM_Score",
           "Angiogenesis_Score", "Phagocytosis_Score",
           "Anti_inflammatory_Score", "Pro_inflammatory_Score")
)

# Inspect metadata
head(scRNA@meta.data)
colnames(scRNA@meta.data)

df <- scRNA@meta.data

# Calculate average scores per macrophage cell type
avg_scores <- df %>%
  group_by(celltype) %>%
  summarise(
    M1_Score = mean(M1_Score1),
    M2_Score = mean(M2_Score2),
    ECM_Score = mean(ECM_Score3),
    Angiogenesis_Score = mean(Angiogenesis_Score4),
    Phagocytosis_Score = mean(Phagocytosis_Score5),
    Anti_inflammatory_Score = mean(Anti_inflammatory_Score6),
    Pro_inflammatory_Score = mean(Pro_inflammatory_Score7)
  ) %>%
  column_to_rownames("celltype")

# Plot heatmap of scaled average scores
pheatmap(
  as.matrix(avg_scores),
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Average Pathway Scores by Celltype"
)

# Save results
write.csv(avg_scores, file = 'Mac_celltype_avg_scores.csv')
