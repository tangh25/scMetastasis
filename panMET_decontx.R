##decontx去除RNA污染 seurat V4
rm(list=ls())
library(Seurat)
library(Seurat)
library(decontX)
# #setwd("E:\\备份\\备份20241021\\2.QC\\0统一基因名")  # 修改为包含 .rds 文件的文件夹路径
# setwd("F:\\coredata\\GEO")
# setwd()
setwd("/home/rstudio/project/Metastasis/data/DISCO/0_rawdata/1_anno/norm/")

#setwd('/home/rstudio/project/Metastasis/data/DISCO/0_rawdata/1_anno/disco_data/norm')
sce<-readRDS("GSE221561_anno_after_qc_norm.rds")
#rownames(sce1)=sce1[,1]
#sce1=sce1[,-1]
# sce <- CreateSeuratObject(counts = sce1, 
#                           min.cells = 3,  
#                           min.features = 200)
# 
# table(Idents(sce))
# check
#p1 <- DimPlot(sce1,label = T)+NoLegend()
#library(devtools)
#devtools::install_github("campbio/decontX")
library(decontX)
library(Seurat)
library(dplyr)
#提取counts矩阵
counts<-sce@assays$RNA@counts
# 得到表达矩阵
#set.seed(123)
#counts <- GetAssayData(object = sce1, layer = "counts")
#counts <- GetAssayData(object = sce1, slot = "counts")
decontX_results<-decontX(counts)
sce$contamination<-decontX_results$contamination
###根据contamination值进行过滤，一般是保留<0.2的细胞,这里过滤掉了7000多个细胞,其实有点多，可能会丢失过多的信息。
sce_filt<-sce[,sce$contamination<0.2]
saveRDS(sce_filt,file="./GSE221561_anno_afterqc_norm_filtered.rds")


########批量处理  seurat v4
library(Seurat)
library(decontX)

# 设置工作目录
#setwd("E:\\备份\\备份20241021\\2.QC\\0统一基因名\\2.GeneID\\norm\\finish")
setwd("F:\\coredata\\GEO")
# 获取所有以_norm.rds结尾的文件
rds_files <- list.files(pattern = "_afterqc_norm.rds$", full.names = TRUE)

# 创建输出文件夹（如果不存在）
output_dir <- "./"
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir)
# }

# 循环处理每个文件
for (rds_file in rds_files) {
  # 读取RDS文件
  sce1 <- readRDS(rds_file)
  #rownames(sce1)=sce1[,1]
  #sce1=sce1[,-1]
  # 创建Seurat对象
  #sce1 <- CreateSeuratObject(counts = sce1, min.cells = 3, min.features = 200)
  
  # 提取counts矩阵
  counts <- sce1@assays$RNA@counts
  
  # 使用decontX去除污染
  decontX_results <- decontX(counts)
  sce1$contamination <- decontX_results$contamination
  
  # 根据contamination值过滤
  sce_filt <- sce1[, sce1$contamination < 0.2]
  
  # 保存过滤后的文件
  output_file <- file.path(output_dir, paste0(basename(rds_file), "_filtered.rds"))
  saveRDS(sce_filt, file = output_file)
  
  # 打印处理完成的消息
  print(paste("Processed and saved:", output_file))
}




###批量操作 Seurat V5
rm(list=ls())
# 加载所需包
library(Seurat)
library(decontX)
library(dplyr)

# 设置工作目录

setwd("/home/rstudio/project/Metastasis/data/GEO/2.QC/0统一基因名/2.GeneID/norm/finish2")  # 修改为包含 .rds 文件的文件夹路径

# 获取所有 .rds 文件的文件名
rds_files <- list.files(pattern = "\\.rds$")

# 遍历每个 .rds 文件
for (file in rds_files) {
  cat("正在处理文件：", file, "\n")
  
  # 读取 RDS 文件
  sce1 <- readRDS(file)
  
  # 确保行名存在
  if (is.null(rownames(sce1))) {
    stop(paste("文件", file, "的行名缺失，请检查数据格式。"))
  }
  
  # 确保数据格式为矩阵
  if (!inherits(sce1, "dgCMatrix")) {
    sce1 <- as(as.matrix(sce1), "dgCMatrix")
  }
  
  # 创建 Seurat 对象
  sce <- CreateSeuratObject(
    counts = sce1, 
    min.cells = 3, 
    min.features = 200
  )
  
  # 提取 counts 矩阵并进行 decontX 处理
  set.seed(123)
  #counts <- GetAssayData(object = sce, layer = "counts")  #seurat V5
  counts <- sce@assays$RNA@counts  #seurat V4
  decontX_results <- decontX(counts)
  
  # 添加污染度信息并过滤细胞
  sce$contamination <- decontX_results$contamination
  sce_filt <- sce[, sce$contamination < 0.2]
  
  # 生成新的文件名
  new_file_name <- sub("\\.rds$", "_rnafilt.rds", file)
  
  # 保存处理后的数据
  saveRDS(sce_filt, file = new_file_name)
  
  cat("已完成：", new_file_name, "\n")
}

cat("所有文件处理完成。\n")


###批量操作 Seurat V4
rm(list=ls())
# 加载所需包
library(Seurat)
library(decontX)
library(dplyr)

# 设置工作目录

setwd("/home/rstudio/project/Metastasis/data/DISCO/0_rawdata/1_anno/GEO25")  # 修改为包含 .rds 文件的文件夹路径

# 获取所有 .rds 文件的文件名
rds_files <- list.files(pattern = "norm\\.rds$")

# 遍历每个 .rds 文件
for (file in rds_files) {
  cat("正在处理文件：", file, "\n")
  
  # 读取 RDS 文件
  sce1 <- readRDS(file)
  
  # 确保行名存在
  if (is.null(rownames(sce1))) {
    stop(paste("文件", file, "的行名缺失，请检查数据格式。"))
  }

  # 创建 Seurat 对象
  sce <- CreateSeuratObject(
    counts = sce1, 
    min.cells = 3, 
    min.features = 200
  )
  
  # 提取 counts 矩阵并进行 decontX 处理
  set.seed(123)
  counts <- sce@assays$RNA@counts  #seurat V4
  decontX_results <- decontX(counts)
  
  # 添加污染度信息并过滤细胞
  sce$contamination <- decontX_results$contamination
  sce_filt <- sce[, sce$contamination < 0.2]
  
  # 生成新的文件名
  new_file_name <- sub("\\.rds$", "_rnafilt.rds", file)
  
  # 保存处理后的数据
  saveRDS(sce_filt, file = new_file_name)
  
  cat("已完成：", new_file_name, "\n")
}

cat("所有文件处理完成。\n")

