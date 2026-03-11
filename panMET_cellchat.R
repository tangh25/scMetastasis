########
## Clear environment and load required libraries
rm(list = ls())
library(Matrix)
library(CellChat)      
library(Seurat)
library(patchwork)
library(BiocParallel)
register(MulticoreParam(workers = 8, progressbar = TRUE))

# ----------------- Load scRNA-seq Seurat object -----------------
scRNA <- readRDS('Pcancer_Fibroblasts_Macrophages_rename.rds')
scRNA$celltype <- scRNA$Subcelltype
Idents(scRNA) <- "celltype" 
table(scRNA$celltype, scRNA$group)

# ----------------- Subset metastatic samples -----------------
cellchat.LS <- subset(scRNA, group == "Metastasis")
scRNA <- cellchat.LS

# ----------------- Prepare normalized expression data -----------------
data.input <- GetAssayData(scRNA, slot = 'data')   
meta <- scRNA@meta.data[, c("orig.ident","celltype")]
colnames(meta) <- c("group","labels")
meta$labels <- gsub(" cells", "", meta$labels)     

# ----------------- Ensure meta and data.input alignment -----------------
identical(rownames(meta), colnames(data.input))

# ----------------- Set cell type order -----------------
meta$labels <- factor(meta$labels, levels = celltype_order)

# Sort meta and data.input according to celltype_order
ordered_indices <- order(meta$labels)
meta <- meta[ordered_indices, ]
data.input <- data.input[, ordered_indices]
identical(rownames(meta), colnames(data.input))

# ----------------- Create CellChat object -----------------
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
levels(cellchat@idents)

# ----------------- Set ligand-receptor database -----------------
CellChatDB <- CellChatDB.human  
CellChatDB.use <- subsetDB(CellChatDB) 
cellchat@DB <- CellChatDB.use

# ----------------- Preprocess CellChat data -----------------
cellchat <- subsetData(cellchat)          
future::plan("multisession", workers = 8) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)  

# ----------------- Infer cell-cell communication -----------------
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net <- subsetCommunication(cellchat)
saveRDS(cellchat,"./cellchat/panMET_MET_Fib_Mac_cellchat_new.rds")
write.csv(df.net,"./cellchat/panMET_MET_Fib_Mac_df_net.csv")

# ----------------- Aggregate and visualize networks -----------------
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

# Circle plots for interactions
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge=FALSE, title.name="Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge=FALSE, title.name="Interaction weights/strength")

# Subset specific cells (Fib ↔ Macro)
fib.cells <- grep("Fib", rownames(cellchat@net$weight), value = TRUE)
macro.cells <- grep("Macro", rownames(cellchat@net$weight), value = TRUE)
mat.fm <- cellchat@net$weight[fib.cells, macro.cells]
mat.mf <- cellchat@net$count[macro.cells, fib.cells]

# Heatmap visualization
pheatmap::pheatmap(mat.mf, border_color = "black", cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 9)

# ----------------- Visualize ligand-receptor contributions -----------------
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR[2,]  # example LR pair

# Individual plots for LR pair
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show,
                     vertex.receiver = vertex.receiver, layout = "hierarchy")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show,
                     layout = "circle")

# ----------------- Outgoing/incoming communication patterns -----------------
selectK(cellchat, pattern = "outgoing")
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 5)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

selectK(cellchat, pattern = "incoming")
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = 6)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")


############# CellChat Analysis Workflow for Two Sample Groups #######################################
rm(list = ls())
library(CellChat)  
library(Matrix)
library(qs)
library(Seurat)
library(patchwork)
library(BiocParallel)
register(MulticoreParam(workers = 8, progressbar = TRUE))

# Load scRNA-seq data
scRNA <- readRDS('Pcancer_Fibroblasts_Macrophages_rename.rds')
scRNA$celltype <- scRNA$Subcelltype
Idents(scRNA) <- "celltype"
table(scRNA$celltype, scRNA$group)

# Split data by sample group
scRNA_left <- scRNA[, scRNA$group %in% c("Metastasis")]
scRNA_right <- scRNA[, scRNA$group %in% c("Primary")]
DimPlot(scRNA_left) | DimPlot(scRNA_right)  # Quick check

####### Create CellChat object for left group
data.input <- GetAssayData(scRNA_left, slot = 'data')  
meta <- scRNA_left@meta.data[, c("orig.ident", "celltype")]
colnames(meta) <- c("group", "labels")
meta$labels <- gsub(" cells", "", meta$labels)
# Sort by cell type
ordered_indices <- order(meta$labels)
meta <- meta[ordered_indices, ]
data.input <- data.input[, ordered_indices]
identical(rownames(meta), colnames(data.input))

# Initialize CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
CellChatDB.use <- subsetDB(CellChatDB.human)  # Use only protein-protein signals
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 8)  # parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net_left <- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat_left <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 

####### Create CellChat object for right group
data.input <- GetAssayData(scRNA_right, slot = 'data')
meta <- scRNA_right@meta.data[, c("orig.ident", "celltype")]
colnames(meta) <- c("group", "labels")
meta$labels <- gsub(" cells", "", meta$labels)
meta$labels <- factor(meta$labels, levels = levels(meta$labels))
ordered_indices <- order(meta$labels)
meta <- meta[ordered_indices, ]
data.input <- data.input[, ordered_indices]
identical(rownames(meta), colnames(data.input))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat@DB <- subsetDB(CellChatDB.human)
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 1)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net_right <- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat_right <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 

### Merge two groups
object.list <- list(left = cellchat_left, right = cellchat_right)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Compare overall interaction number and strength
compareInteractions(cellchat, show.legend = F, group = c(1, 2))
compareInteractions(cellchat, show.legend = F, group = c(1, 2), measure = "weight")

# Compare interactions between cell types
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat, measure = "weight")

# Circle plot of interaction counts
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T,
                   label.edge = F, edge.weight.max = weight.max[2],
                   edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# Scatter plot of signaling roles
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- lapply(1:length(object.list), function(i) {
  netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
})
patchwork::wrap_plots(plots = gg)

# Differential signaling from left to right
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:7), comparison = c(1, 2))
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:7), comparison = c(1, 2),
                        max.dataset = 2, title.name = "Increased signaling in Left", angle.x = 45)
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:7), comparison = c(1, 2),
                        max.dataset = 1, title.name = "Decreased signaling in Right", angle.x = 45)
gg1 + gg2

# Identify up/down-regulated ligand-receptor pairs
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = "right",
                                      features.name = "right.merged", only.pos = FALSE, thresh.pc = 0.1,
                                      thresh.fc = 0.05, thresh.p = 0.05)
net <- netMappingDEG(cellchat, features.name = "right.merged", variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "left", ligand.logFC = 0.05)
net.down <- subsetCommunication(cellchat, net = net, datasets = "right", ligand.logFC = -0.05)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# Bubble plots for up/down signaling
pairLR.use.up <- net.up[, "interaction_name", drop = F]
pairLR.use.down <- net.down[, "interaction_name", drop = F]
netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11),
                 comparison = c(1,2), angle.x = 90, remove.isolate = T)
netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11),
                 comparison = c(1,2), angle.x = 90, remove.isolate = T)

# Violin plot for gene expression
cellchat@meta$datasets <- factor(cellchat@meta$datasets, levels = c("left", "right"))
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T, type = "violin")











