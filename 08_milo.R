rm(list = ls())
############################
# Load required packages
############################
library(stringr)
library(Seurat)
library(miloR)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(qs)
library(BiocParallel)

# Parallel computation
register(MulticoreParam(workers = 20, progressbar = TRUE))
############################
# Load annotated dataset
load('data_RunUMAP.RData')
scRNA <- sc_data
rm(sc_data)

############################
# Check metadata
############################

table(scRNA$orig.ident)
table(scRNA$anno_Sub)

############################
# Convert Seurat object to SingleCellExperiment
############################

scRNA_pre <- as.SingleCellExperiment(scRNA)

reducedDimNames(scRNA_pre)

############################
# Construct Milo object
############################

scRNA_milo <- Milo(scRNA_pre)

reducedDim(scRNA_milo,"UMAP") <- reducedDim(scRNA_pre,"UMAP")

rm(scRNA_pre)
rm(scRNA)

scRNA_milo

############################
# Build KNN graph
############################

scRNA_milo <- buildGraph(scRNA_milo,
                         k = 30,
                         d = 15,
                         reduced.dim = "PCA")

############################
# Define representative neighborhoods
############################

scRNA_milo <- makeNhoods(scRNA_milo,
                         prop = 0.2,
                         k = 30,
                         d = 15,
                         refined = TRUE,
                         reduced_dims = "PCA")

plotNhoodSizeHist(scRNA_milo)

############################
# Count cells per neighborhood
############################

scRNA_milo <- countCells(scRNA_milo,
                         meta.data = as.data.frame(colData(scRNA_milo)),
                         sample = "sample")

head(nhoodCounts(scRNA_milo))

############################
# Experimental design matrix
############################

traj_design <- data.frame(colData(scRNA_milo))[,c("sample", "group")]

traj_design <- distinct(traj_design)

rownames(traj_design) <- traj_design$sample

traj_design <- traj_design[colnames(nhoodCounts(scRNA_milo)), , drop=FALSE]

traj_design$group <- factor(traj_design$group,
                            level = c("Metastasis", "Primary"))

traj_design

############################
# Calculate neighborhood distances
############################

scRNA_milo <- calcNhoodDistance(scRNA_milo,
                                d = 15,
                                reduced.dim = "PCA")

saveRDS(scRNA_milo,file = 'Fib_milo2.rds')

############################
# Differential abundance testing
############################

da_results <- testNhoods(scRNA_milo,
                         design = ~ group,
                         design.df = traj_design)

da_results %>%
  arrange(SpatialFDR) %>%
  head()

write.csv(da_results,file='da_results.csv')

############################
# Basic visualization
############################

ggplot(da_results, aes(PValue)) +
  geom_histogram(bins = 50)

ggplot(da_results,
       aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = 1)

############################
# Build neighborhood graph
############################

scRNA_milo <- buildNhoodGraph(scRNA_milo)

############################
# UMAP visualization
############################

umap_pl <- plotReducedDim(scRNA_milo,
                          dimred = "UMAP",
                          colour_by="group",
                          text_by = "anno_Sub",
                          text_size = 3,
                          point_size=0.5) +
  guides(fill="none")

############################
# Plot differential abundance neighborhoods
############################

nh_graph_pl <- plotNhoodGraphDA(scRNA_milo,
                                da_results,
                                layout="UMAP",
                                alpha=0.9)

umap_pl + nh_graph_pl +
  plot_layout(guides="collect")

############################
# Annotate neighborhoods by cell type
############################

da_results <- annotateNhoods(scRNA_milo,
                             da_results,
                             coldata_col = "anno_Sub")

head(da_results)

plotDAbeeswarm(da_results,
               group.by = "anno_Sub",
               alpha = 0.9)

write.csv(da_results,file = 'da_results.csv')

############################
# Identify DA neighborhood groups
############################

scRNA_milo <- buildNhoodGraph(scRNA_milo)

da_results <- groupNhoods(scRNA_milo,
                          da_results,
                          max.lfc.delta = 2,
                          da.fdr = 0.9)

plotNhoodGroups(scRNA_milo,
                da_results,
                layout="UMAP")

plotDAbeeswarm(da_results,
               "NhoodGroup",
               alpha = 0.9)

############################
# Identify marker genes
############################

keep.rows <- rowSums(logcounts(scRNA_milo)) != 0
scRNA_milo <- scRNA_milo[keep.rows, ]

dec <- modelGeneVar(scRNA_milo)

hvgs <- getTopHVGs(dec, n = 2000)

############################
# Differential expression between neighborhood groups
############################

nhood_markers <- findNhoodGroupMarkers(
  scRNA_milo,
  da_results,
  subset.row = hvgs,
  aggregate.samples = TRUE,
  sample_col = "orig.ident"
)

head(nhood_markers)

############################
# Visualization of differential genes
############################

ggplot(nhood_markers,
       aes(logFC_2,-log10(adj.P.Val_2))) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_hline(yintercept = 0.05)

markers <- nhood_markers$GeneID[
  nhood_markers$adj.P.Val_2 < 0.05 &
  nhood_markers$logFC_2 > 1
]

plotNhoodExpressionGroups(
  scRNA_milo,
  da_results,
  features = markers,
  subset.nhoods = da_results$NhoodGroup %in% c('2','5'),
  scale = TRUE,
  grid.space = "fixed"
)

############################
# Differential expression within specific neighborhoods
############################

dge_5 <- testDiffExp(
  scRNA_milo,
  da_results,
  design = ~ location,
  meta.data = data.frame(colData(scRNA_milo)),
  subset.row = rownames(scRNA_milo),
  subset.nhoods = da_results$NhoodGroup == "5"
)

dge_5
