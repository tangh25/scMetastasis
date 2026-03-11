# Load necessary libraries
library(nichenetr)        # NicheNet analysis
library(Seurat)           # Single-cell RNA-seq analysis
library(SeuratObject)
library(tidyverse)        # Data manipulation and visualization

# Load Seurat object and update to current version
seuratObj <- readRDS("seuratObj.rds")
seuratObj <- UpdateSeuratObject(seuratObj)
seuratObj <- alias_to_symbol_seurat(seuratObj, "human")  # Convert gene aliases to symbols

# Visualize cells in UMAP/tSNE
DimPlot(seuratObj)

# Set organism for ligand-receptor networks
organism <- "human"
if(organism == "human"){
  lr_network <- readRDS("lr_network_human_21122021.rds")
  ligand_target_matrix <- readRDS("ligand_target_matrix_nsga2r_final.rds")
  weighted_networks <- readRDS("weighted_networks_nsga2r_final.rds")
} else if(organism == "mouse"){
  lr_network <- readRDS("./mouse/lr_network_mouse_21122021.rds")
  ligand_target_matrix <- readRDS("./mouse/ligand_target_matrix_nsga2r_final_mouse.rds")
  weighted_networks <- readRDS("./mouse/weighted_networks_nsga2r_final_mouse.rds")
}

# Remove duplicate ligand-receptor pairs
lr_network <- lr_network %>% distinct(from, to)

# Define receiver cell type and get expressed genes
receiver <- "Macro_LYVE1"
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)

# Identify expressed receptors and potential ligands for receiver
all_receptors <- unique(lr_network$to)
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

# Define sender cell types
sender_celltypes <- c("Fib_HOPX")

# Get expressed genes in sender cells
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

# Focus ligands that are expressed in both sender and receiver
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Define conditions for differential expression analysis
condition_oi <- "Metastasis"
condition_reference <- "Primary"

# Subset receiver cells and compute DE genes
seurat_obj_receiver <- subset(seuratObj, idents = receiver)
DE_table_receiver <- FindMarkers(
  object = seurat_obj_receiver,
  ident.1 = condition_oi,
  ident.2 = condition_reference,
  group.by = "aggregate",
  min.pct = 0.05
) %>% rownames_to_column("gene")

# Select DE genes as target gene set
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

# Define background genes expressed in receiver
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

# Predict ligand activities using NicheNet
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

# Rank ligands by activity (AUPR)
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))

# Plot histogram of ligand activities, highlight Top20 cutoff
p_hist_lig_activity <- ggplot(ligand_activities, aes(x = aupr_corrected)) +
  geom_histogram(color = "black", fill = "darkorange") +
  geom_vline(aes(xintercept = min(ligand_activities %>% top_n(20, aupr_corrected) %>% pull(aupr_corrected))),
             color = "red", linetype = "dashed", size = 1) +
  labs(x = "Ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

# Select Top20 ligands based on AUPR
best_upstream_ligands <- ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

# Heatmap of ligand activities for Top20 ligands
LigandActivity <- function(ligand_activities, best_upstream_ligands){
  vis_ligand_aupr <- ligand_activities %>% 
    filter(test_ligand %in% best_upstream_ligands) %>%
    column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% as.matrix()
  
  heatmapplot <- make_heatmap_ggplot(vis_ligand_aupr, "Prioritized ligands", "Ligand activity", legend_title = "AUPR", color = "darkorange") +
    theme(axis.text.x.top = element_blank())
  
  return(list(heatmapplot = heatmapplot, vis_ligand_aupr_ncol = ncol(vis_ligand_aupr)))
}
LigandActivity(ligand_activities, best_upstream_ligands)

# Prepare ligand-target regulatory links
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(active_ligand_target_links_df, ligand_target_matrix, cutoff = 0.33)

# Heatmap of ligand-target regulatory potential
RegulatoryPlot <- function(best_upstream_ligands, active_ligand_target_links){
  order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets <- intersect(unique(active_ligand_target_links_df$target), rownames(active_ligand_target_links))
  vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])
  
  regheatmap <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                    color = "purple", legend_title = "Regulatory potential") +
    scale_fill_gradient2(low = "whitesmoke", high = "purple")
  return(list(regheatmap = regheatmap, vis_ligand_targetncol = ncol(vis_ligand_target)))
}
RegulatoryPlot(best_upstream_ligands, active_ligand_target_links)

# Infer ligand-receptor interactions for Top ligands
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(best_upstream_ligands, expressed_receptors, lr_network, weighted_networks$lr_sig)

# Heatmap of ligand-receptor interaction potential
receptorplot <- function(ligand_receptor_links_df, best_upstream_ligands){
  vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(ligand_receptor_links_df, best_upstream_ligands, order_hclust = "both")
  make_heatmap_ggplot(t(vis_ligand_receptor_network), y_name = "Ligands", x_name = "Receptors",
                      color = "mediumvioletred", legend_title = "Prior interaction potential")
}
receptorplot(ligand_receptor_links_df, best_upstream_ligands)

# Sender-focused analysis: filter ligands expressed in sender cells
ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(20, aupr_corrected) %>% pull(test_ligand) %>% unique()

# Heatmap of sender-focused ligand activity
p_ligand_auprs <- LigandActivity(ligand_activities, best_upstream_ligands)
p_ligand_aupr <- p_ligand_auprs$heatmapplot
p_ligand_aupr

# Heatmap of ligand-target regulatory links (sender-focused)
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 100) %>%
  bind_rows() %>% drop_na()
active_ligand_target_links <- prepare_ligand_target_visualization(active_ligand_target_links_df, ligand_target_matrix, cutoff = 0.33)
p_ligand_targets <- RegulatoryPlot(best_upstream_ligands, active_ligand_target_links)
p_ligand_target <- p_ligand_targets$regheatmap
p_ligand_target

# Heatmap of ligand-receptor interactions (sender-focused)
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(best_upstream_ligands, expressed_receptors, lr_network, weighted_networks$lr_sig)
p_ligand_receptor <- receptorplot(ligand_receptor_links_df, best_upstream_ligands)
p_ligand_receptor

# DotPlot of Top ligand expression in sender cells
p_dotplot <- DotPlot(subset(seuratObj, celltype %in% sender_celltypes), features = rev(best_upstream_ligands), cols = "RdYlBu") +
  coord_flip() + scale_y_discrete(position = "right")
p_dotplot

# Heatmap of log fold-change of Top ligands in sender cells
celltype_order <- levels(Idents(seuratObj))
DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype,
  seurat_obj = seuratObj,
  condition_colname = "aggregate",
  condition_oi = condition_oi,
  condition_reference = condition_reference,
  celltype_col = "celltype",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands
) %>% reduce(full_join) %>% column_to_rownames("gene")
vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])
p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc, "Prioritized ligands", "LFC in Sender",
                                        low_color = "midnightblue", mid_color = "white",
                                        mid = median(vis_ligand_lfc), high_color = "red",
                                        legend_title = "LFC")
p_lfc

# Compare sender-agnostic vs sender-focused ligand activities
make_line_plot(ligand_activities = ligand_activities_all, potential_ligands = potential_ligands_focused) +
  theme(plot.title = element_text(size = 11, hjust = 0.1, margin = margin(0,0,-5,0)))

# Combine all NicheNet plots into a single figure
figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 12), axis.text.y = element_text(size = 9),
                    axis.text.x = element_text(size = 9, angle = 90, hjust = 0)) + ylab("Expression in Sender"),
  p_lfc + theme(legend.position = "none", axis.title.y = element_blank()),
  p_ligand_target + theme(legend.position = "none", axis.title.y = element_blank()),
  align = "hv",
  nrow = 1,
  rel_widths = c((p_ligand_auprs$vis_ligand_aupr_ncol)+6, ncol(vis_ligand_lfc)+7, ncol(vis_ligand_lfc)+8, (p_ligand_targets$vis_ligand_targetncol)+8)
)
legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
  nrow = 1,
  align = "h",
  rel_widths = c(1.5,1,1,1)
)
combined_plot <- cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
ggsave("allplot.png", width = 28, height = 16, plot = combined_plot)
combined_plot
