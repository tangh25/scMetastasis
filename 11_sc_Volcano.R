## After cell annotation, visualize DEGs for each cell subtype
library(Seurat)
library(jjVolcano)
library(dplyr)
library(scRNAtoolVis)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)  # human gene annotation

setwd("./data")
# ----------------- Load Seurat object -----------------
sce <- readRDS('./data_umap.rds')
sce <- snhx  # optional overwrite
ncol(sce)  # number of cells
ncol(sce@assays$RNA@counts)   

# ----------------- Visualize variable genes -----------------
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave("Variablegenes.pdf", plot = plot2, width = 8, height = 5, dpi = 300)

# ----------------- UMAP plots -----------------
DimPlot(sce, reduction = "umap", group.by = "anno_sub")  # colored by cell subtype
DimPlot(sce, reduction = "umap", group.by = "anno_sub", split.by = "group")  # split by group

Idents(sce) <- "anno_sub"

# ----------------- Find marker genes -----------------
markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "markers_jjVolcano.csv")

# ----------------- Rename columns to match Python DEG outputs -----------------
markers <- read.csv('diff.csv')
colnames(markers)[colnames(markers) == "logfoldchanges"] <- "avg_log2FC"
colnames(markers)[colnames(markers) == "pvals"] <- "p_val"
colnames(markers)[colnames(markers) == "pvals_adj"] <- "p_val_adj"
markers$cluster <- markers$anno_Major

# ----------------- Filter significant DEGs -----------------
markers2 <- markers %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
write_csv(markers2, file='sig_diff.csv')
head(markers2)

# ----------------- Volcano plot -----------------
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(20)  # generate 20 distinct colors

jjVolcano(diffData = markers2, tile.col = mycolors)  # basic volcano
jjVolcano(diffData = markers2, topGeneN = 5, col.type = 'adjustP', tile.col = mycolors)  # top 5 genes

# Polar (circular) volcano plot
jjVolcano(diffData = markers, tile.col = mycolors, polar = TRUE)

# ----------------- Visualize cell proportion per sample -----------------
cellRatio <- cellRatioPlot(object = sce,
                           sample.name = "orig.ident",
                           celltype.name = "celltype",
                           flow.curve = 0.5,
                           fill.col = mycolors) +
             theme(axis.text.x = element_text(angle = 45, hjust = 1))
cellRatio

# ----------------- KEGG enrichment for upregulated genes -----------------
gid <- bitr(unique(markers$gene), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
DEG <- full_join(markers, gid, by = c('gene' = 'SYMBOL'))

KEGG <- compareCluster(ENTREZID ~ cluster, data = DEG, fun = 'enrichKEGG')
dotplot(KEGG, label_format = 40) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradient(high = "#4b5cc4", low = "#FE8D3C")

# ----------------- Save Seurat object -----------------
save(sce, file = 'scRNA_deg.rdata')
