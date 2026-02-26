# ==============================
# MISTy spatial interaction analysis
# ==============================

# remotes::install_github("saezlab/mistyR")

library(mistyR)     # MISTy spatial modeling framework
library(Seurat)     # Spatial Seurat object handling
library(distances)  # Distance calculations
library(future)     # Parallel computing support
library(tidyverse)  # Data manipulation
library(recipes)

# Load spatial Seurat object
GBM4 <- readRDS('adata_vis_GSM5420753_region.rds')

# Visualize spatial distribution
SpatialPlot(GBM4)

########################################################
# Extract cell-type composition (e.g., deconvolution output)
########################################################

GBM4@meta.data

# Columns to exclude (technical or metadata fields)
exclude_cols <- c("region", "coord_key", "nCount_SCT", "nFeature_SCT",
                  'SCT_snn_res.0.5','sample','_indices','_scvi_batch',
                  "orig.ident", "nCount_Spatial", "nFeature_Spatial",
                  'array_row','array_col','_scvi_labels',
                  "_index", "in_tissue",'seurat_clusters')

# Keep only composition-related columns
composition <- GBM4@meta.data[, !(colnames(GBM4@meta.data) %in% exclude_cols)]

head(composition)

# Make column names syntactically valid (required by MISTy)
name_map <- data.frame(
  original = colnames(composition),
  safe = make.names(colnames(composition), unique = TRUE)
)

colnames(composition) <- name_map$safe

########################################################
# Extract spatial coordinates
########################################################

geometry <- GetTissueCoordinates(
  GBM4,
  cols = c("imagerow", "imagecol"),
  scale = NULL
)

########################################################
# Define MISTy views
########################################################

# Compute distance matrix between spots
geom_dist <- as.matrix(distances(geometry))

# Distance to nearest neighbor for each spot
dist_nn <- apply(geom_dist, 1, function(x) (sort(x)[2]))

# Define paraview radius (mean NN distance + SD)
paraview_radius <- ceiling(mean(dist_nn + sd(dist_nn)))

# Create intrinsic (intraview) and spatial (paraview)
GBM_views <- create_initial_view(composition) %>%
  add_paraview(geometry, l = paraview_radius, family = "gaussian")

########################################################
# Run MISTy model
########################################################

library(ridge)

run_misty(
  GBM_views,
  "st_data"
)

########################################################
# Downstream analysis
########################################################

# Collect MISTy results
misty_results <- collect_results(
  "st_data"
)

# Plot variance explained improvement by multi-view modeling
misty_results %>%
  plot_improvement_stats("multi.R2") %>%
  plot_improvement_stats("gain.R2")

########################################################
# Interaction importance visualization
########################################################

# Heatmap of intrinsic (intraview) interactions
misty_results %>%
  plot_interaction_heatmap(view = "intra", clean = FALSE)

# Identify predictors for specific cell types
misty_results$importances.aggregated %>%
  filter(view == "intra", Predictor == "Fib_HOPX") %>%
  arrange(-Importance) %>%
  print(n = Inf)

# Filter interactions using gain.R2 threshold
misty_results %>%
  plot_interaction_heatmap(
    view = "intra",
    clean = FALSE,
    trim = 0.05,
    trim.measure = "gain.R2",
    cutoff = 0.5
  )

# Save results
saveRDS(misty_results, file = './GSM5420753_misty_results.rds')

# Reload if needed
misty_results <- readRDS('./GSM5420753_misty_results.rds')

########################################################
# Visualize spatial distribution of selected cell types
########################################################

SpatialFeaturePlot(
  GBM4,
  features = c("Fib_HOPX", "Macro_LYVE1"),
  image.alpha = 1,
  pt.size.factor = 1.6
)

SpatialFeaturePlot(
  GBM4,
  features = c("Fib_HOPX", "Macro_LYVE1"),
  image.alpha = 0,
  pt.size.factor = 1.6
)

# Save spatial object
save(GBM4, file = './GSM5420753_misty.rdata')

load('./GSM5420753_misty.rdata')

########################################################
# Co-localization analysis (custom visualization)
########################################################

library(dplyr)
library(ggplot2)

# Extract interaction importance for a specific predictor
tcell_interactions <- misty_results$importances.aggregated %>%
  filter(view == "intra", Predictor == "Fib_HOPX") %>%
  arrange(-Importance)

# Cell types to analyze
celltypes_to_analyze <- c("Fib_HOPX", "Macro_LYVE1")

########################################################
# Function to plot spatial co-localization
########################################################

plot_colocalization <- function(celltype1, celltype2,
                                composition_df,
                                geometry_df,
                                interaction_df) {
  
  # Co-occurrence score (product of proportions)
  cooc <- composition_df[, celltype1] *
          composition_df[, celltype2]
  
  # Spatial dataframe
  spatial_df <- data.frame(
    x = geometry_df$imagecol,
    y = geometry_df$imagerow,
    cooc = cooc
  )
  
  # Extract interaction importance
  interaction_imp <- interaction_df %>%
    filter(Target == celltype2) %>%
    pull(Importance)
  
  # Plot spatial co-localization
  p <- ggplot(spatial_df,
              aes(x = x, y = y, color = cooc)) +
    geom_point(size = 2) +
    scale_color_gradientn(
      colors = c('blue', 'yellow', 'red'),
      name = 'CO-loc_Score'
    ) +
    scale_y_reverse() +  # Match spatial orientation
    labs(
      title = paste(celltype1, "&", celltype2,
                    "co-localization"),
      subtitle = ifelse(
        length(interaction_imp) > 0,
        paste("Interaction Importance:",
              round(interaction_imp[1], 3)),
        "Interaction not calculated or found"
      )
    ) +
    theme_classic() +
    coord_fixed()
  
  return(p)
}

########################################################
# Generate co-localization plots
########################################################

co_loc_plots <- list()

for (ct in setdiff(celltypes_to_analyze, "HOPX..Fib")) {
  
  co_loc_plots[[ct]] <- plot_colocalization(
    celltype1 = "HOPX..Fib",
    celltype2 = ct,
    composition_df = composition,
    geometry_df = geometry,
    interaction_df = tcell_interactions
  )
}

# Display plots
co_loc_plots
