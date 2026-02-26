#BiocManager::install("genefilter")
#BiocManager::install("genefilter")

#devtools::install_github("JEFworks/MUDAN")
#devtools::install_github("data2intelligence/SpaCET")
Sys.setenv(RETICULATE_PYTHON="/home/htang/.conda/envs/TumorBoundary/bin/python.exe")
library(SpaCET)
#空转数据的读取
visiumPath = file.path(system.file(package = "SpaCET"), "extdata/Visium_BC")

SpaCET_obj = create.SpaCET.object.10X(visiumPath = visiumPath)
#SpaCET的内部结构如下，input存储着输入数据信息，其中counts存储着表达矩阵的原始counts数，spotCoordinates存储着spot的空间坐标，image存储着图像信息。results部分存储分析结果，目前还没有分析结果。
SpaCET_obj@input$counts[1:8,1:6]
#质控
SpaCET_obj = SpaCET.quality.control(SpaCET_obj)
#可视化UMI counts数和基因在空间中的分布
SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "QualityControl",
  spatialFeatures = c("UMI", "Gene"),
  imageBg = T
)

#反卷积
SpaCET_obj = SpaCET.deconvolution(SpaCET_obj, cancerType = "BRCA", coreNo = 8)
SpaCET_obj@results$deconvolution$propMat[1:13, 1:6]

#可视化
SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "CellFraction",
  spatialFeatures = c("Malignant", "Macrophage")
)

#一次性可视化所有一级和二级细胞类型的空间分布
SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "CellFraction",
  spatialFeatures = "All",
  sameScaleForFraction = T,
  pointSize = 0.1,
  nrow = 5
)

#交互式空间可视化功能
SpaCET.visualize.spatialFeature(SpaCET_obj,interactive=TRUE)


#细胞通讯
#细胞共定位分析
SpaCET_obj = SpaCET.CCI.colocalization(SpaCET_obj)

SpaCET.visualize.colocalization(SpaCET_obj)

#计算L-R互作网络
SpaCET_obj = SpaCET.CCI.LRNetworkScore(SpaCET_obj, coreNo = 8)

SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "LRNetworkScore",
  spatialFeatures = c("Network_Score", "Network_Score_pv")
)

SpaCET_obj = SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair = c("CAF", "Macrophage M2"))
SpaCET.visualize.cellTypePair(SpaCET_obj, cellTypePair = c("CAF", "Macrophage M2"))

#肿瘤组织边界分析
# Identify the Tumor-Stroma Interface
SpaCET_obj = SpaCET.identify.interface(SpaCET_obj)

# Visualize the Interface
SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", spatialFeatures = "Interface")

#可视化细胞对共定位和肿瘤组织边界
# Combine the interface and interaction spots
SpaCET_obj = SpaCET.combine.interface(SpaCET_obj, cellTypePair = c("CAF", "Macrophage M2"))

# Visualize the Interface. 
SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", spatialFeatures = "Interface&CAF_Macrophage M2")

SpaCET.distance.to.interface(SpaCET_obj, cellTypePair = c("CAF", "Macrophage M2"))

#肿瘤异质性亚群的揭示
SpaCET_obj = SpaCET.deconvolution.malignant(SpaCET_obj, coreNo = 8)
SpaCET_obj@results$deconvolution$propMat[c("Malignant cell state A", "Malignant cell state B"), 1:6]

SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "CellFraction",
  spatialFeatures = c("Malignant", "Malignant cell state A", "Malignant cell state B"),
  nrow = 1
)

write_rds(SpaCET_obj, "brac.deconv.rds")


###########
#单细胞数据辅助反卷积


PDAC_Path = system.file("extdata", "oldST_PDAC", package = 'SpaCET')
load(paste0(PDAC_Path, "/st_PDAC.rda"))

counts[1:6,1:5]
spotCoordinates[1:5, ]


SpaCET_obj = create.SpaCET.object(
  counts,
  spotCoordinates,
  imagePath = NA,    #imagePath指定的是H&E染色的图像，这套数据没有，所以是NA。
  platform = "oldST"    #platform是空间组学平台，可以是"Visium"、"OldST"或者"Slide-seq"。"OldST"是早期原位空间捕获技术
)


# calculate the QC metrics
SpaCET_obj = SpaCET.quality.control(SpaCET_obj)

# plot the QC metrics
SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "QualityControl",
  spatialFeatures = c("UMI", "Gene"),
  imageBg = T
)

#单细胞数据的导入
#单细胞转录组数据需要首先处理成sc_counts和sc_annotation的形式，这两个矩阵就是大家非常熟悉的counts矩阵和metadata，很容易从seurat对象得到
sc_counts[1:6, 1:5]

sc_annotation[1:6,]

#这个软件非常好用的地方还在于，它可以同时对大群和亚群进行反卷积，只需要把大群对应的亚群存储在sc_lineageTree这个列表里就行了。
head(sc_lineageTree)

#反卷积
SpaCET_obj = SpaCET.deconvolution.matched.scRNAseq(
  SpaCET_obj,
  sc_counts = sc_counts,
  sc_annotation = sc_annotation,
  sc_lineageTree = sc_lineageTree,
  coreNo = 8
)

SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "CellFraction",
  spatialFeatures = c("Cancer clone A", "Cancer clone B", "Acinar cells", "Ductal - CRISP3 high/centroacinar like"),
  nrow = 2
)

SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "GeneExpression",
  spatialFeatures = c("TM4SF1", "S100A4", "PRSS1", "CRISP3"),
  nrow = 2
)

write_rds(SpaCET_obj, "brac.deconv.sc.rds")

########基因集评分
#类似于GSVA和AUCell，SpaCET也能计算基因集评分，这里我就省略前面数据导入和质控的过程，直接计算几个基因集的分数。
#Hallmark,Hallmark是Msigdb数据库中的一个经典基因集
SpaCET_obj = SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "Hallmark")
SpaCET_obj@results$GeneSetScore[1:6, 1:6]

SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "GeneSetScore",
  spatialFeatures = c("HALLMARK_HYPOXIA", "HALLMARK_TGF_BETA_SIGNALING")
)

#Cancer cell state
#文献：https://www.nature.com/articles/s41588-022-01141-9 Cancer cell state基因集记录的是肿瘤细胞的各种状态对应的基因，可以采用相同的方式进行基因集评分和可视化。
# run gene set calculation
SpaCET_obj = SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "CancerCellState")

# visualize two gene sets
SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "GeneSetScore",
  spatialFeatures = c("CancerCellState_Cycle", "CancerCellState_cEMT")
)
#TLS
#文献：https://www.researchsquare.com/article/rs-3921508 TLS就是一个单独的基因集，它包含有30个三级淋巴结构的标志基因，计算这个基因集的评分并进行可视化，可以明显看出三级淋巴结构的位置。
SpaCET_obj = SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "TLS")

SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "GeneSetScore",
  spatialFeatures = "TLS"
)
#自定义基因集
gmt1 = list(
  Tcell = c("CD2", "CD3E", "CD3D"),
  Myeloid = c("SPI1", "FCER1G", "CSF1R")
)
SpaCET_obj = SpaCET.GeneSetScore(SpaCET_obj, GeneSets = gmt1)
SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "GeneSetScore",
  spatialFeatures = "Tcell"
)

#在Msigdb官网上能够下载到各种通路的gmt文件，而gmt文件读取进来就是列表格式，所以能直接进行分析。这里举的例子是读取KEGG通路进行评分分析。
gmt2 = read.gmt("~/ref/human/msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
SpaCET_obj = SpaCET.GeneSetScore(SpaCET_obj, GeneSets = gmt2)
SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "GeneSetScore",
  spatialFeatures = "KEGG_DNA_REPLICATION"
)




