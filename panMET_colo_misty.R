#MISTy
### 实战5:空间共定位分析
# remotes::install_github("saezlab/mistyR")
library(mistyR)
#library(CARD)
library(Seurat)
library(mistyR)
library(distances)
library(future)
library(tidyverse) 
library(recipes)

#setwd("/home/rstudio/project/Metastasis/data/ST/LUAD/GSE179572/GSM5420753")
setwd("/home/rstudio/project/Metastasis/data/ST/LUAD/GSE179572/GSM5420753")

# 载入实战2保存的 GBM4 空转数据
#load('GBM4.rdata')
GBM4 <- readRDS('adata_vis_GSM5420753_region.rds')
SpatialPlot(GBM4)


# 载入实战3保存的 CARD反卷积数据
#load('CARD_obj.rdata')

# # 提取反卷积细胞成分
GBM4@meta.data
# 指定需要排除的列名
exclude_cols <- c("region", "coord_key", "nCount_SCT", "nFeature_SCT",'SCT_snn_res.0.5','sample','_indices','_scvi_batch',
                  "orig.ident", "nCount_Spatial", "nFeature_Spatial",'array_row','array_col','_scvi_labels',
                  "_index", "in_tissue",'seurat_clusters')

# 提取剩余的列作为 composition
composition <- GBM4@meta.data[, !(colnames(GBM4@meta.data) %in% exclude_cols)]

# 查看前几行
head(composition)
#统一处理列名
name_map <- data.frame(
  original = colnames(composition),
  safe = make.names(colnames(composition), unique = TRUE)
)

colnames(composition) <- name_map$safe

# composition <- as.data.frame(CARD_obj@Proportion_CARD)
# #将列名(细胞名)中的空格替换为下划线,否则后续会报错
# colnames(composition) <- gsub(" ", "_", colnames(composition), fixed = TRUE)
# 提取空间位置信息
geometry <- GetTissueCoordinates(GBM4, cols = c("imagerow", "imagecol"), scale = NULL)

## 首先，需要定义一个intraview，以捕捉一个点内的细胞类型比例，
## 为了捕捉周围组织中细胞类型比例的分布，我们添加了一个paraview，
## 我们选择的半径是到最近邻居的距离的平均值加上标准偏差，
## 使用family=gaussian 计算每个点的权重，然后运行MISTy并收集结果。

#Calculating the radius
geom_dist <- as.matrix(distances(geometry))  
dist_nn <- apply(geom_dist, 1, function(x) (sort(x)[2]))
paraview_radius <- ceiling(mean(dist_nn+ sd(dist_nn)))

# Create views
GBM_views <- create_initial_view(composition) %>%
  add_paraview(geometry, l= paraview_radius, family = "gaussian")

# Run misty and collect results
library(ridge)
run_misty(GBM_views, "/home/rstudio/project/Metastasis/data/ST/LUAD/GSE179572/GSM5420753")

##下游分析,读取上一步运行结果
misty_results <- collect_results("/home/rstudio/project/Metastasis/data/ST/LUAD/GSE179572/GSM5420753")
#与intraview相比，周围组织的细胞类型在多大程度上可以解释斑点的细胞类型组成？在这里，我们可以看到两种不同的统计数据：multi.R2 显示了由多视图模型解释的总方差；gain.R2显示了从全景到全景的可解释方差的增加。
misty_results %>%
  plot_improvement_stats("multi.R2") %>% 
  plot_improvement_stats("gain.R2")

#为了解释些贡献的具体关系，我们可以分别可视化每种细胞类型在预测每个视图的细胞类型分布中的重要性。

##首先, 对于intrinsic view
## (view = "intra", clean = F)
misty_results %>% plot_interaction_heatmap(view = "intra", clean = F)


##计算哪些细胞可以预测OPC_like细胞
misty_results$importances.aggregated %>%
  filter(view == "intra", Predictor == "HOPX..Fib") %>%
  arrange(-Importance) %>%
  print(n = Inf)   # 显示所有行

misty_results$importances.aggregated %>%
  filter(view == "intra", Predictor == "Macro_C1QC") %>%
  arrange(-Importance) %>%
  print(n = Inf)   # 显示所有行


##指定paraview_radius，设置闸值(trim 默认 1.75%) Paraview 在很大程度上解释了其分配，Oligo细胞的最佳预测指标是OPC_like细胞
unique(misty_results$importances.aggregated$view)

paraview_radius
misty_results %>% plot_interaction_heatmap(view = "intra", clean = F, 
                                           trim = 0.05, trim.measure = "gain.R2",
                                           cutoff = 0.5)
saveRDS(misty_results,file = './GSM5420753_misty_results.rds')
misty_results <- readRDS('./GSM5420753_misty_results.rds')
## 反卷积结果加至ST对象
#GBM4[['CARD']] <- CreateAssayObject(counts = t(composition))
##可视化细胞类型空间分布
SpatialFeaturePlot(GBM4, features = c("HOPX+ Fib","Macro_LYVE1"), image.alpha = 1,pt.size.factor = 1.6)   #image.alpha = 0或1
SpatialFeaturePlot(GBM4, features = c("HOPX+ Fib","Macro_LYVE1"), image.alpha = 0,pt.size.factor = 1.6)   #image.alpha = 0或1
save(GBM4,file = './GSM5420753_misty.rdata')
load('./GSM5420753_misty.rdata')

##############提取特定细胞
#' library(dplyr)
#' library(ggplot2)
#' # library(ggpubr) # 如果需要使用 ggarrange 来排列图表
#' 
#' # --- 1. 修正的初始化和交互性计算 ---
#' 
#' # 修正：确保使用正确的管道符(%>%)，并使用半角逗号(,)分隔 filter 内部条件。
#' tcell_interactions <- misty_results$importances.aggregated %>%
#'   filter(view == "intra", Predictor=="HOPX..Fib") %>% 
#'   arrange(-Importance)
#' 
#' # 修正细胞类型列表的变量名拼写
#' celltypes_to_analyze <- c("HOPX..Fib", "Macro_LYVE1", "Macro_C1QC") 
#' 
#' # --- 2. 修正后的共定位绘图函数 (plot_colocalization) ---
#' 
#' #' 绘制两种细胞类型的空间共定位图
#' #'
#' #' @param celltype1 预测细胞类型 (Predictor) 名称
#' #' @param celltype2 目标细胞类型 (Target) 名称
#' #' @param composition_df 包含细胞类型分数或比例的矩阵/数据框
#' #' @param geometry_df 包含空间坐标 (imagecol, imagerow) 的数据框
#' #' @param interaction_df 包含 MISTy 交互作用重要性得分的数据框
#' plot_colocalization <- function(celltype1, celltype2, composition_df, geometry_df, interaction_df){
#'   
#'   # 1. 计算共定位得分 (Co-occurrence Score)
#'   # 修正了 composition[,celltypel] 的拼写错误为 composition_df[, celltype1]
#'   cooc <- composition_df[, celltype1] * composition_df[, celltype2]
#'   
#'   # 2. 获取空间信息 (修正了数据框创建时的逗号和 geometry$imagerow 的引用)
#'   spatial_df <- data.frame(
#'     x = geometry_df$imagecol,
#'     y = geometry_df$imagerow,
#'     cooc = cooc
#'   )
#'   
#'   # 3. 获取相互作用重要性值 (修正了管道符和 Target 筛选逻辑)
#'   interaction_imp <- interaction_df %>%
#'     filter(Target == celltype2) %>% 
#'     pull(Importance)
#'   
#'   # 4. 绘图 (修正了颜色拼写、标题拼写和标点符号)
#'   p <- ggplot(spatial_df, aes(x = x, y = y, color = cooc)) +
#'     geom_point(size = 2) +
#'     scale_color_gradientn(
#'       colors = c('blue', 'yellow', 'red'),
#'       name = 'CO-loc_Score' # 修正了未闭合的括号
#'     ) +
#'     labs(
#'       title = paste(celltype1, "&", celltype2, "co-localization"),
#'       # 修正了副标题中的中文全角逗号
#'       subtitle = ifelse(length(interaction_imp) > 0,
#'                         paste("Interaction Importance:", round(interaction_imp[1], 3)),
#'                         "Interaction not calculated or found")
#'     ) +
#'     theme_classic() +
#'     coord_fixed() # 保持坐标轴比例
#'   
#'   return(p)
#' }
#' 
#' # --- 3. 生成共定位图的循环示例 ---
#' 
#' co_loc_plots <- list()
#' # 循环遍历 'celltypes_to_analyze' 中除了 'HOPX+ Fib' 之外的所有细胞类型
#' for (ct in setdiff(celltypes_to_analyze, "HOPX..Fib")){
#'   # 确保将实际的数据对象传递给函数：composition, geometry, tcell_interactions
#'   co_loc_plots[[ct]] <- plot_colocalization(
#'     celltype1 = "HOPX..Fib",
#'     celltype2 = ct, # 目标细胞类型：Macro_LYVE1 和 Macro_C1QC
#'     composition_df = composition,
#'     geometry_df = geometry,
#'     interaction_df = tcell_interactions
#'   )
#' }
#' 
#' # 使用 ggpubr::ggarrange 或 patchwork 库来排列图表 (如果需要)
#' #ggarrange(plotlist = co_loc_plots, ncol = 2, nrow = 1) 
#' co_loc_plots
#' # 保存 MISTy 结果 (保持原意)
#' save(misty_results, composition, geometry, file = 'LUAD/GSE179572/GSM5420753_misty_full_results.rdata')
#' #load('./LUAD/GSE179572/GSM5420753_misty_full_results.rdata')
#' # load('/home/rstudio/project/Metastasis/data/ST/LUAD/GSE179572/GSM5420753/LUAD/GSE179572/GSM5420753_misty.rdata')

###########################翻转y轴
library(dplyr)
library(ggplot2)

# --- 1. 提取特定细胞的交互重要性 ---
tcell_interactions <- misty_results$importances.aggregated %>%
  filter(view == "intra", Predictor == "HOPX..Fib") %>% 
  arrange(-Importance)

# 指定要分析的细胞类型
celltypes_to_analyze <- c("HOPX..Fib", "Macro_LYVE1", "Macro_C1QC") 

# --- 2. 共定位绘图函数 ---
plot_colocalization <- function(celltype1, celltype2, composition_df, geometry_df, interaction_df){
  
  # 计算共定位得分 (Co-occurrence Score)
  cooc <- composition_df[, celltype1] * composition_df[, celltype2]
  
  # 获取空间信息
  spatial_df <- data.frame(
    x = geometry_df$imagecol,
    y = geometry_df$imagerow,
    cooc = cooc
  )
  
  # 获取相互作用重要性值
  interaction_imp <- interaction_df %>%
    filter(Target == celltype2) %>% 
    pull(Importance)
  
  # 绘图
  p <- ggplot(spatial_df, aes(x = x, y = y, color = cooc)) +
    geom_point(size = 2) +
    scale_color_gradientn(
      colors = c('blue', 'yellow', 'red'),
      name = 'CO-loc_Score'
    ) +
    scale_y_reverse() +  # 翻转 y 轴，使图像方向与原始空间图一致
    labs(
      title = paste(celltype1, "&", celltype2, "co-localization"),
      subtitle = ifelse(length(interaction_imp) > 0,
                        paste("Interaction Importance:", round(interaction_imp[1], 3)),
                        "Interaction not calculated or found")
    ) +
    theme_classic() +
    coord_fixed() # 保持 x、y 坐标比例
  
  return(p)
}

# --- 3. 循环生成共定位图 ---
co_loc_plots <- list()
for (ct in setdiff(celltypes_to_analyze, "HOPX..Fib")){
  co_loc_plots[[ct]] <- plot_colocalization(
    celltype1 = "HOPX..Fib",
    celltype2 = ct,
    composition_df = composition,
    geometry_df = geometry,
    interaction_df = tcell_interactions
  )
}

# --- 4. 可选：使用 ggpubr 或 patchwork 排列图表 ---
# library(ggpubr)
# ggarrange(plotlist = co_loc_plots, ncol = 2, nrow = 1)

co_loc_plots

# 保存 MISTy 结果 (保持原意)
save(misty_results, composition, geometry, file = './GSM5420753_misty_full_results.rdata')
#load('./GSM5420753_misty_full_results.rdata')
# load('/home/rstudio/project/Metastasis/data/ST/LUAD/GSE179572/GSM5420753/LUAD/GSE179572/GSM5420753_misty.rdata')

# geometry_df$imagecol → x 坐标
# 
# 表示水平位置（column），从左到右增加
# 
# geometry_df$imagerow → y 坐标
# 
# 表示垂直位置（row），从上到下增加（在原始图像中）
# 
# 所以 x 和 y 就是每个 spot 或 cell 在 原始空间组织切片图像上的位置。

# 为什么需要翻转 y
# 
# 在 R 的 ggplot2 中：
# 
# (0,0) 默认在左下角
# 
# y 值向上增加
# 
# 而在大多数空间转录组原始图像中（比如 10x Visium 或 imager 数据）：
# 
# (0,0) 在左上角
# 
# y 值向下增加
# 
# 所以直接用 y = imagerow 会导致图像 上下翻转。
# 解决方法就是在 ggplot 中加 scale_y_reverse()，或者在数据层面做 y = max(imagerow) - imagerow。