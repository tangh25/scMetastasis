########
rm(list = ls())
##http://172.25.51.2:8716/
rm(list = ls())
library(Matrix)
library(CellChat) # 2.1.2版本
#library(qs)
library(Seurat)
library(patchwork)
library(BiocParallel)
register(MulticoreParam(workers = 8, progressbar = TRUE))
scRNA <- readRDS('Pcancer_Fibroblasts_Macrophages_rename.rds')
scRNA$celltype <- scRNA$Subcelltype
Idents(scRNA) <- "celltype" 
#rm(data_harmony)
table(scRNA$celltype,scRNA$group)

# 获取 Metastasis 
cellchat.LS <- subset(scRNA, group == "Metastasis")

scRNA <- cellchat.LS
unique(scRNA$group)
#scRNA <- cellchat.NL
##CellChat要求输入标准化后的表达数据----
data.input <- GetAssayData(scRNA, assay = "RNA", layer = "data")
table(scRNA@meta.data$celltype)
Idents(scRNA) <- scRNA$celltype
DimPlot(scRNA,label = F)


#######创建Cellchat对象
data.input <- GetAssayData(scRNA, slot = 'data') # normalized data matrix
meta <- scRNA@meta.data[,c("orig.ident","celltype")]
colnames(meta) <-  c("group","labels")
table(meta$labels)
meta$labels <- gsub(" cells", "", meta$labels)
#meta$labels <- sub("\\(.*\\)", "", meta$labels)
table(meta$labels)

identical(rownames(meta),colnames(data.input))
## 根据研究情况进行细胞排序
celltype_order <- c(
  "Fib_ADAMDEC1 (lamina propria)", 
  "Fib_ADH1B (alveolar)", 
  "Fib_CD74 (ap)", 
  "Fib_COL15A1", 
  "Fib_CTNNB1", 
  "Fib_HGF", 
  "Fib_HHIP(SMC)", 
  "Fib_HOPX", 
  "Fib_HSPA6",
  "Fib_IL6", 
  "Fib_LRRC15", 
  "Fib_Mesothelial", 
  "Fib_MMP1", 
  "Fib_MYH11(SMC)", 
  "Fib_Pericyte", 
  "Fib_PI16", 
  "Fib_PRG4 (synovial lining)", 
  "Fib_SOX6 (epithelial crypt)",
  "Fib_SPRP2", 
  "Fib_STMN1 (proliferative)",
  "Macro_C1QC", 
  "Macro_FN1",
  "Macro_GPNMB", 
  "Macro_IL1B",
  "Macro_INHBA", 
  "Macro_ISG15",
  "Macro_LYVE1", 
  "Macro_NLRP3",
  "Macro_SPP1"
)
meta$labels <- factor(meta$labels ,levels = celltype_order)

table(meta$labels)

# 根据 meta$labels 的顺序进行排序
ordered_indices <- order(meta$labels)
# 对 meta 和 data.input 进行排序
meta <- meta[ordered_indices, ]
data.input <- data.input[, ordered_indices]
identical(rownames(meta),colnames(data.input))

# 构建cellchat
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
levels(cellchat@idents)
#  [1] "T/NK"              "Th1"               "Th17"              "Tm"                "Treg"             
#  [6] "Naive T"           "ELK4+T"            "ZNF793+T"          "ZSCAN12+T"         "B"                
# [11] "VSMCs"             "endothelial"       "epithelial/cancer" "fibroblasts"       "mast"             
# [16] "myeloid"           "plasma"            "proliferative"   


#####设置配体-受体相互作用数据库
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

# 除“非蛋白信号”外，使用所有CellChatDB数据进行细胞-细胞通信分析
CellChatDB.use <- subsetDB(CellChatDB)

# 使用所有CellChatDB数据进行细胞-细胞通信分析
# 研究者不建议以这种方式使用它，因为CellChatDB v2包含“非蛋白信号”(即代谢和突触信号)。
# CellChatDB.use <- CellChatDB 

# 在构建的cellchat中设定需要使用的数据库
cellchat@DB <- CellChatDB.use

###预处理细胞-细胞通讯分析的表达数据
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole databa
future::plan("multisession", workers = 1) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#默认情况下,cellchat使用object@data.signaling进行网络推断
#同时也提供了projectData函数,通过扩散过程基于高置信度实验验证的蛋白质互作网络中的邻近节点对基因表达值进行平滑处理。该功能在处理测序深度较浅的单细胞数据时尤为有用，因其能减少信号基因（特别是配体/受体亚基可能存在的零表达）的dropout效应。不担心其可能在扩散过程引入伪影，因其仅会引发极微弱的通讯信号。
# 原来是projectData，新版是smoothData函数
cellchat <- smoothData(cellchat, adj = PPI.human)

####细胞-细胞通信网络的推理
# 该分析的关键参数是类型，即计算每个细胞组的平均基因表达的方法。默认情况下，type = “triMean”，产生较少但更强的交互。当设置 type = “truncatedMean” 时，应为trim分配一个值，从而产生更多交互。请详细检查上述计算每个细胞组平均基因表达的方法。
# 使用的是投射到PPI网络的模式时候需要用FALSE。如果使用了raw data就需要设置为TRUE
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE) 

# 如果所研究的信号没有被测到，可以采用如下函数进行探查，trim设为0.1或者0.05
# computeAveExpr(cellchat, features = c("CXCL12","CXCR4"),type =  "truncatedMean",trim = 0.1)
# 如果发现修改参数之后所研究的信号被测到了，那就修改代码如下
# cellchat <- computeCommunProb(cellchat, type =  "truncatedMean",trim = 0.1,raw.use = FALSE) 

# min.cells是设置阈值，最小是需要10个细胞参与通讯推断(可以自定义)
cellchat <- filterCommunication(cellchat, min.cells = 10)

#在信号通路水平上推断细胞间通讯
# CellChat通过汇总与每个信号通路相关的所有配体-受体相互作用的通信概率来计算信号通路水平上的通信概率。 
# NB:推断的每个配体-受体对的细胞间通信网络和每个信号通路分别存储在槽'net'和'netP'中。
cellchat <- computeCommunProbPathway(cellchat)

#数据提取，subsetCommunication函数，一般全部提取并保存
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #表示从细胞群 1 和 2 向细胞群 4 和 5 推断出的细胞间通讯。
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df.net <- subsetCommunication(cellchat)

saveRDS(cellchat,"./cellchat/panMET_MET_Fib_Mac_cellchat_new.rds")
cellchat <- readRDS("../cellchat/panMET_MET_Fib_Mac_cellchat_new.rds")
save(df.net,file = "./cellchat/panMET_MET_Fib_Mac_df_net.Rdata")
write.csv(df.net,"./cellchat/panMET_MET_Fib_Mac_df_net.csv")

###可视化流程
#CellChat可以可视化聚合的蜂窝间通信网络。circle展示互作数目
# 计算聚合细胞-细胞通信网络
# 互作网络整合,可以设置soure和target，不设置就是默认全部
cellchat <- aggregateNet(cellchat)
# 可视化 这两张图分别展示了互作的数量和互作的权重。其中每个颜色代表了不同的细胞，箭头代表了顺序，线的粗细代表了数量/权重。
groupSize <- as.numeric(table(cellchat@idents)) 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


###circle plot 指定 source / target
#先子集化矩阵，再画图
cells.use <- c(
  "Fib_HOPX",
  "Macro_C1QC", 
  "Macro_FN1",
  "Macro_GPNMB", 
  "Macro_IL1B",
  "Macro_INHBA", 
  "Macro_ISG15",
  "Macro_LYVE1", 
  "Macro_NLRP3",
  "Macro_SPP1"
)

mat <- cellchat@net$weight[cells.use, cells.use]
netVisual_circle(
  mat,
  sources.use = "Fib_HOPX",
  targets.use = setdiff(cells.use, "Fib_HOPX"),
  vertex.weight = rowSums(mat),
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Fib_HOPX → Macrophage interactions"
)

mat2 <- cellchat@net$count[cells.use, cells.use]
netVisual_circle(
  mat2,
  sources.use = "Fib_HOPX",
  targets.use = setdiff(cells.use, "Fib_HOPX"),
  vertex.weight = rowSums(mat),
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Fib_HOPX → Macrophage interactions"
)

netVisual_circle(
  mat2,
  sources.use = "Fib_HOPX",
  targets.use = setdiff(cells.use, "Fib_HOPX"),
  vertex.weight = rowSums(mat2),
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of LR interactions"
)


##heatmap展示互作数据
pheatmap::pheatmap(cellchat@net$count, border_color = "black", 
                   cluster_cols = F, fontsize = 10, cluster_rows = F,
                   display_numbers = F,number_color="black",number_format = "%.0f")

##只展示 Fib ↔ Macro 之间的互作，必须构造一个“行 = Fib，列 = Macro（或反之）的矩阵”，而不是 Fib+Macro 的方阵
fib.cells <- c(
  "Fib_ADAMDEC1 (lamina propria)", 
  "Fib_ADH1B (alveolar)", 
  "Fib_CD74 (ap)", 
  "Fib_COL15A1", 
  "Fib_CTNNB1", 
  "Fib_HGF", 
  "Fib_HHIP(SMC)", 
  "Fib_HOPX", 
  "Fib_HSPA6",
  "Fib_IL6", 
  "Fib_LRRC15", 
  "Fib_Mesothelial", 
  "Fib_MMP1", 
  "Fib_MYH11(SMC)", 
  "Fib_Pericyte", 
  "Fib_PI16", 
  "Fib_PRG4 (synovial lining)", 
  "Fib_SOX6 (epithelial crypt)",
  "Fib_SPRP2", 
  "Fib_STMN1 (proliferative)"
)

macro.cells <- c(
  "Macro_C1QC", 
  "Macro_FN1",
  "Macro_GPNMB", 
  "Macro_IL1B",
  "Macro_INHBA", 
  "Macro_ISG15",
  "Macro_LYVE1", 
  "Macro_NLRP3",
  "Macro_SPP1"
)

mat.fm <- cellchat@net$weight[fib.cells, macro.cells]

mat.mf <- cellchat@net$count[macro.cells, fib.cells]

pheatmap::pheatmap(
  mat.mf,
  border_color = "black",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize = 9
)






##贝克图
#展示每一个celltype作为source与其他celltype的互作情况
#指定顺序和指定颜色
celltype_order <- c(
  "Fib_ADAMDEC1 (lamina propria)", 
  "Fib_ADH1B (alveolar)", 
  "Fib_CD74 (ap)", 
  "Fib_COL15A1", 
  "Fib_CTNNB1", 
  "Fib_HGF", 
  "Fib_HHIP(SMC)", 
  "Fib_HOPX", 
  "Fib_HSPA6",
  "Fib_IL6", 
  "Fib_LRRC15", 
  "Fib_Mesothelial", 
  "Fib_MMP1", 
  "Fib_MYH11(SMC)", 
  "Fib_Pericyte", 
  "Fib_PI16", 
  "Fib_PRG4 (synovial lining)", 
  "Fib_SOX6 (epithelial crypt)",
  "Fib_SPRP2", 
  "Fib_STMN1 (proliferative)",
  "Macro_C1QC", 
  "Macro_FN1",
  "Macro_GPNMB", 
  "Macro_IL1B",
  "Macro_INHBA", 
  "Macro_ISG15",
  "Macro_LYVE1", 
  "Macro_NLRP3",
  "Macro_SPP1"
)


# 将颜色向量命名为矩阵的行名

mat <- as.data.frame(cellchat@net$weight)
mat <- mat[celltype_order,]#行排序
mat <- mat[,celltype_order] %>% as.matrix()
# 生成颜色向量（例如使用彩虹色）
names(color.use) <- rownames(mat)
color.use <- rainbow(nrow(mat))

# 1. 检查您的细胞类型名称 (CellChat通常使用这个)
cell_type_names <- rownames(mat)

# 2. 检查您的颜色向量长度是否与细胞类型数量匹配
# 您提供的颜色数量是 30 个
provided_colors <- c(
  "#FF0000", "#FF3500", "#FF6A00", "#FF9E00", "#FFD300", 
  "#F6FF00", "#C1FF00", "#8DFF00", "#58FF00", "#23FF00", 
  "#00FF12", "#00FF46", "#00FF7B", "#00FFB0", "#00FFE5", 
  "#00E5FF", "#00B0FF", "#007BFF", "#0046FF", "#0012FF", 
  "#2300FF", "#5800FF", "#8D00FF", "#C100FF", "#F600FF", 
  "#FF00D3", "#FF009E", "#FF006A", "#FF0035" # 这里一共是 29 个颜色
) 

# 如果您的细胞类型数量大于 29，provided_colors 会不足。
if (length(cell_type_names) > length(provided_colors)) {
  stop("提供的颜色数量不足以分配给所有的细胞类型。")
}

# 3. 创建命名向量 (这是解决错误的关键步骤)
# 我们只取与细胞类型数量相匹配的颜色
color.use.named <- setNames(
  provided_colors[1:length(cell_type_names)], 
  cell_type_names
)

# 4. 运行绘图循环 (使用修正后的命名向量)
par(mfrow = c(5,4), xpd=TRUE, mar = c(1, 1, 1, 1))

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2, 
    vertex.weight = groupSize, 
    weight.scale = T, 
    arrow.size = 0.05,
    arrow.width = 1, 
    edge.weight.max = max(mat), 
    title.name = rownames(mat)[i],
    color.use = color.use.named # <<-- 使用命名向量
  )
}

# # 如果图片显示不全,需要考虑是不是重新设置mfrow   其中每个颜色代表了不同的细胞，箭头代表了顺序，线的粗细代表了数量/权重。
# par(mfrow = c(5,4), xpd=TRUE,mar = c(1, 1, 1, 1))
# for (i in 1:nrow(mat)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, 
#                    weight.scale = T, arrow.size=0.05,
#                    arrow.width=1, edge.weight.max = max(mat), 
#                    title.name = rownames(mat)[i],
#                    color.use = color.use)
# }

##层次结构图
#层次图：vertex.receiver设定了source的细胞类型(实心圆圈)，空心圆圈代表target。左半边图片先把代表vertex.receiver的圆圈放在了中间，显示了不同细胞类型对其他细胞的作用情况，右半边图片把代表其他细胞的圆圈放在了中间，显示了不同细胞类型对其余细胞的作用情况，线的粗细代表互作的强度。
cellchat@netP$pathways
levels(cellchat@idents) 
pathways.show <- "COLLAGEN"  #COLLAGEN,CXCL 

# Hierarchy plot
# vertex.receiver定义层次图的左边细胞
vertex.receiver = seq(1:9) # a numeric vector
netVisual_aggregate(cellchat, signaling = pathways.show,
                    vertex.receiver = vertex.receiver,layout= "hierarchy")
# vertex.size = groupSize)  
# circle plot
#Circle图和分组弦图：不同颜色代表不同细胞类型，结合箭头代表作用方向，粗细代表互作的强度，两个图很表达方式几乎一样。
netVisual_aggregate(cellchat, signaling = pathways.show,layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# 分组弦图
#分组弦图：笔者把T细胞相关的设定成一组之后会跟其他细胞之间存在一定的分组关系
group.cellType <- c(rep("T/NK", 9), "B","VSMCs","endothelial","epithelial/cancer",
                    "fibroblasts","mast","myeloid","plasma","proliferative" )
levels(cellchat@idents)
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show,
                     group = group.cellType,
                     title.name = paste0(pathways.show, " signaling network"))

# heatmap
# 热图：
# 
# 行： 热图的行代表了不同的细胞类型，这些细胞作为信号的发送者。
# 列： 热图的列代表了不同的细胞类型，这些细胞作为信号的接收者。
# 颜色深浅： 热图中的颜色深浅表示了通讯概率的大小。颜色越深，表示通讯概率越高，这意味着发送方细胞和接收方细胞之间的信号传递越强。
# 通信概率（Communication Prob.）： 右侧的颜色条是颜色映射的参考。图中的深红色表示较高的通讯概率（靠近0.05），浅色表示较低的通讯概率（靠近0或更低）。
# 顶部数值范围 (0 - 0.3)： 显示不同细胞类型接受到总传入信号强度(incoming）
# 右侧部分的数值范围 (0 - 0.4)：显示不同细胞类型发出的总传出信号强度(outcoming）

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")


##计算配-受体对信号通路的贡献并可视化
# 计算配-受体对信号通路的贡献并可视化
netAnalysis_contribution(cellchat, signaling = pathways.show)

# 可视化由单个配体-受体对介导的细胞间通讯
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show,
                            geneLR.return = FALSE)
pairLR
#   interaction_name
# 1      CXCL1_CXCR2
# 2      CXCL2_CXCR2
# 3      CXCL3_CXCR2
# 4      CXCL8_CXCR2
# 5     CXCL12_CXCR4
# 6     CXCL12_ACKR3
# 7     CXCL16_CXCR6
LR.show <- pairLR[2,] # show one ligand-receptor pair

# Hierarchy plot
vertex.receiver = seq(1,9) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  
                     pairLR.use = LR.show, 
                     vertex.receiver = vertex.receiver,
                     layout = "hierarchy")

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, 
                     pairLR.use = LR.show, layout = "circle")

##基于配-受体结果进一步可视化
#Bubble plot, 使用netVisual_bubble显示一些单元组与其他单元组之间的所有重要相互作用

# 需要指定source和target
# sources.use是发出信号的细胞系,target.use是接受信号的细胞系
levels(cellchat@idents) 
netVisual_bubble(cellchat, sources.use = seq(21:29), 
                 targets.use = c(8), remove.isolate = FALSE)
ggsave("bubbleplot_nont.pdf",width = 7,height = 20)


netVisual_bubble(cellchat, sources.use = c(8), 
                 targets.use = c(27), remove.isolate = FALSE)
ggsave("bubbleplot_nont.pdf",width = 7,height = 20)
# 还可以增加signaling参数用于展示特定的配受体   !!!!!!!!!!!!!!!!!!
cellchat@netP$pathways
netVisual_bubble(cellchat, sources.use = c(8), 
                 targets.use = c(21:29), 
                 signaling = c("COLLAGEN"),
                 remove.isolate = FALSE)

ggsave("bubbleplot2.pdf",width = 2,height = 10)




# 自定义signaling输入展示-所有通路汇总之后
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CXCL","COLLAGEN","MK","CD99",
                                                        "LAMININ","APP","DESMOSOME"))
netVisual_bubble(cellchat, sources.use = c(8),
                 targets.use = c(21:29), 
                 pairLR.use = pairLR.use,
                 remove.isolate = TRUE)
ggsave("bubbleplot-LR.pdf",width = 5,height = 15)

# 可以通过增加下面的参数去设置X轴上的顺序
# sort.by.target = T
# sort.by.source = T
# sort.by.source = T, sort.by.target = T
# sort.by.source = T, sort.by.target = T, sort.by.source.priority = FALSE

####Chord diagram
# 这里进行绘制时建议设定singaling或者减少互作的细胞
# 内容大多会导致不出图
cellchat@netP$pathways
netVisual_chord_gene(cellchat, sources.use = c(1:9), 
                     targets.use = c(13), 
                     signaling = c("VEGF"),
                     lab.cex = 0.5,
                     legend.pos.y = 30)

# 显示从某些细胞组(sources.use)到其他细胞组(targets.use)的所有重要信号通路。
netVisual_chord_gene(cellchat, sources.use = c(1:9), 
                     targets.use = c(13), 
                     slot.name = "netP", 
                     legend.pos.x = 10)
# NB: Please ignore the note when generating the plot such as “Note: The first link end is drawn out of sector ‘MIF’.”. If the gene names are overlapped, you can adjust the argument small.gap by decreasing the value.

##使用小提琴/点图绘制信号转导基因表达分布

# CellChat可以使用Seurat包装函数plotGeneExpression绘制与L-R对或信号通路相关的信号转导基因的基因表达分布。
# 该功能提供 “violin”、“dot”、“bar” 三种类型的可视化。
# 或用户可以使用 extractEnrichedLR 提取与推断的 L-R 对或信号通路相关的信号转导基因，然后使用Seurat或其他软件包绘制基因表达。
plotGeneExpression(cellchat, signaling = "VEGF", 
                   enriched.only = TRUE, 
                   type = "violin")

##计算并可视化网络中心性得分
cellchat@netP$pathways
pathways.show <- "VEGF"
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, 
                                          slot.name = "netP") 
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                  width = 8, height = 2.5, font.size = 10)

##行（Sender, Receiver, Mediator, Influencer）：行表示在信号通路网络中，不同细胞类型扮演的角色：
# 
# Sender：信号的发送者，即哪些细胞类型是主要的信号发出者。
# Receiver：信号的接收者，即哪些细胞类型是主要的信号接收者。
# Mediator：中介者，用于识别在信号传播路径中起中介作用的细胞群体。
# Influencer：影响者，用于识别在整个网络中对信息传播影响最大的细胞群体。 列（不同的细胞类型）：列表示具体的细胞类型。颜色深浅（Importance）：颜色的深浅表示每个细胞类型在特定角色中的重要性。颜色越深，表示该细胞类型在这个角色中的重要性越高（例如信号传递的强度或频率越大）；颜色越浅，表示该细胞类型在这个角色中的重要性较低。
# 

##在二维空间中可视化占优势的发送者(源)和接收者(目标)
# 从所有信号通路对聚合细胞-细胞通信网络的信号作用分析
gg1 <- netAnalysis_signalingRole_scatter(cellchat);gg1
# 对特定细胞间通讯网络的信号作用分析
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL"));gg2
gg1 + gg2

##识别对某些细胞群的输出或输入信号贡献最大的信号
# ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
# ht1
# ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
# ht2
# ht1 + ht2
# class(ht1)

# 特定的signaling
cellchat@netP$pathways
htout <- netAnalysis_signalingRole_heatmap(cellchat, 
                                           pattern = "outgoing",
                                           signaling = c("ICAM","TGFb"))
htout

htcome <- netAnalysis_signalingRole_heatmap(cellchat, 
                                            pattern = "incoming",
                                            signaling = c("ICAM","TGFb"))
htcome

######识别整体通信模式/以探索多种细胞类型和信号通路如何协调运作
#outgoing
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")

# 当输出模式的数量为5时，Cophenetic值和Silhouette值都开始突然下降。
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "outgoing", 
                                          k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")



selectK(cellchat, pattern = "incoming")

# 当输出模式的数量为6时，Cophenetic值和Silhouette值都开始突然下降。
nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "incoming", 
                                          k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")





#############两组样本Cellchat分析流程#######################################
rm(list = ls())
library(CellChat) # 2.1.2版本
library(Matrix)
library(qs)
library(Seurat)
library(patchwork)
library(BiocParallel)
register(MulticoreParam(workers = 8, progressbar = TRUE))
scRNA <- readRDS('Pcancer_Fibroblasts_Macrophages_rename.rds')
scRNA$celltype <- scRNA$Subcelltype
Idents(scRNA) <- "celltype" 
#rm(data_harmony)
table(scRNA$celltype,scRNA$group)

# 拆分数据
scRNA_left <- scRNA[,scRNA$group%in% c("Metastasis")]
scRNA_right <-scRNA[,scRNA$group%in% c("Primary")]
# check
DimPlot(scRNA_left)|DimPlot(scRNA_right)

#######分别创建Cellchat对象
##cellchat_left
# 自己分析的时候记得改名，这里是scRNA_left
data.input <- GetAssayData(scRNA_left, slot = 'data') # normalized data matrix
meta <- scRNA_left@meta.data[,c("orig.ident","celltype")]
colnames(meta) <-  c("group","labels")
table(meta$labels)
meta$labels <- gsub(" cells", "", meta$labels)
#meta$labels <- sub("\\(.*\\)", "", meta$labels)
table(meta$labels)

identical(rownames(meta),colnames(data.input))
## 根据研究情况进行细胞排序
celltype_order <- c(
  "T/NK", 
  "Th1", 
  "Th17", 
  "Tm", 
  "Treg", 
  "Naive T", 
  "ELK4+T", 
  "ZNF793+T", 
  "ZSCAN12+T",
  "B", 
  "VSMCs", 
  "endothelial", 
  "epithelial/cancer", 
  "fibroblasts", 
  "mast", 
  "myeloid", 
  "plasma", 
  "proliferative"
)
meta$labels <- factor(meta$labels ,levels = celltype_order)
table(meta$labels)

# 根据 meta$labels 的顺序进行排序
ordered_indices <- order(meta$labels)
# 对 meta 和 data.input 进行排序
meta <- meta[ordered_indices, ]
data.input <- data.input[, ordered_indices]
identical(rownames(meta),colnames(data.input))

# 构建cellchat
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
CellChatDB <- CellChatDB.human  
#除“非蛋白信号”外，使用所有CellChatDB数据进行细胞-细胞通信分析
CellChatDB.use <- subsetDB(CellChatDB)
# 使用所有CellChatDB数据进行细胞-细胞通信分析,研究者不建议以这种方式使用它，因为CellChatDB v2包含“非蛋白信号”(即代谢和突触信号)。
#CellChatDB.use <- CellChatDB 

# 在构建的cellchat中设定需要使用的数据库
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole databa
future::plan("multisession", workers = 1) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net_left <- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat_left <- netAnalysis_computeCentrality(cellchat, 
                                               slot.name = "netP") 

##cellchat_right
# 自己分析的时候记得改名，这里是scRNA_left
data.input <- GetAssayData(scRNA_right, slot = 'data') # normalized data matrix
meta <- scRNA_right@meta.data[,c("orig.ident","celltype")]
colnames(meta) <-  c("group","labels")
table(meta$labels)
meta$labels <- gsub(" cells", "", meta$labels)
#meta$labels <- sub("\\(.*\\)", "", meta$labels)
table(meta$labels)

identical(rownames(meta),colnames(data.input))
## 根据研究情况进行细胞排序
celltype_order <- c(
  "T/NK", 
  "Th1", 
  "Th17", 
  "Tm", 
  "Treg", 
  "Naive T", 
  "ELK4+T", 
  "ZNF793+T", 
  "ZSCAN12+T",
  "B", 
  "VSMCs", 
  "endothelial", 
  "epithelial/cancer", 
  "fibroblasts", 
  "mast", 
  "myeloid", 
  "plasma", 
  "proliferative"
)
meta$labels <- factor(meta$labels ,levels = celltype_order)
table(meta$labels)

# 根据 meta$labels 的顺序进行排序
ordered_indices <- order(meta$labels)
# 对 meta 和 data.input 进行排序
meta <- meta[ordered_indices, ]
data.input <- data.input[, ordered_indices]
identical(rownames(meta),colnames(data.input))

# 构建cellchat
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
# human
cellchatDB <- CellChatDB.human  
#除“非蛋白信号”外，使用所有CellChatDB数据进行细胞-细胞通信分析
CellChatDB.use <- subsetDB(CellChatDB)
# 使用所有CellChatDB数据进行细胞-细胞通信分析,研究者不建议以这种方式使用它，因为CellChatDB v2包含“非蛋白信号”(即代谢和突触信号)。
#CellChatDB.use <- CellChatDB 

# 在构建的cellchat中设定需要使用的数据库
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole databa
future::plan("multisession", workers = 1) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net_right <- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat_right <- netAnalysis_computeCentrality(cellchat, 
                                                slot.name = "netP") 

###整合两组样本
object.list <- list(left=cellchat_left,right=cellchat_right)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

##比较相互作用的总数和互作强度
##整体比较
gg1 <- compareInteractions(cellchat, show.legend = F, 
                           group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, 
                           group = c(1,2),measure = "weight")
gg1 + gg2

##比较不同细胞群体之间的互作数量和强度
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

##显示两个数据集中不同细胞群体间相互作用数量或相互作用强度的热图
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

##Circle图显示多个数据集中不同细胞群体之间的相互作用数量或相互作用强度
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

##Circle图显示粗细胞类型之间的相互作用数量或相互作用强度差异
# 按照celltype_order进行分组
celltype_order
# 分组，这里笔者瞎分一下，没有任何生物学意义哈
group.cellType <- c(rep("T/NK", 9), rep("Group1", 4), rep("Group2", 5))
group.cellType <- factor(group.cellType, levels = c("T/NK", "Group1", "Group2"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
# interactions
par(mfrow = c(1,2), xpd=TRUE)
for (i in1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# weight.merged
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)


# 在二维空间中比较主要来源和目标
# 识别发送或接收信号发生显著变化的细胞群体
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

patchwork::wrap_plots(plots = gg)


### 可视化从left样本到right样本的发出与接收信号差异变化
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Th1", signaling.exclude = "MIF")

gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Treg", signaling.exclude = c("MIF"))

patchwork::wrap_plots(plots = list(gg1,gg2))

##比较每个信号通路或配体-受体对的总体信息流。
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

##比较与每个细胞群体相关的传出和传入信号模式
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

###识别上调和下调的信号配体-受体对
#通过比较通信概率来识别功能障碍的信号
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:7),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:7),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Left", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:7),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in right", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2


# 通过差异表达分析识别功能障碍性信号
# 对每个细胞群体在两种生物条件之间进行差异表达分析，然后根据sender细胞中配体的倍数变化和接收细胞中受体的倍数变化来获得上调和下调的信号。

# 设定"实验组"
pos.dataset = "right"
features.name = paste0(pos.dataset, ".merged")
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 

net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "left",ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "right",ligand.logFC = -0.05, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# 自定义特征和感兴趣的细胞群体找到所有显著的outgoing/incoming/both向信号
df <- findEnrichedSignaling(object.list[[2]], features = c("CCL19", "CXCL12"), idents = c("Treg", "Tm"), pattern ="outgoing")


###可视化已识别的上调/下调信号配体-受体对
# 气泡图
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

# 和弦图
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

# 词云图
library(wordcloud)
computeEnrichmentScore(net.down, species = 'human', variable.both = TRUE)
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)

##可视化，使用层次图、圆形图或弦图直观比较细胞间通讯
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
# Chord diagram 的另外一种形式
# group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
# names(group.cellType) <- levels(object.list[[1]]@idents)
# pathways.show <- c("CXCL") 
# par(mfrow = c(1,2), xpd=TRUE)
# for (i in 1:length(object.list)) {
#   netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
# }

par(mfrow = c(1, 2), xpd=TRUE)
for (i in1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 4, targets.use = c(5:6), lab.cex = 0.5, title.name = paste0("Signaling from Tm - ", names(object.list)[i]))
}

# compare all the interactions sending from fibroblast to inflamatory immune cells
# par(mfrow = c(1, 2), xpd=TRUE)
# for (i in 1:length(object.list)) {
#   netVisual_chord_gene(object.list[[i]], sources.use = c(1,2, 3, 4), targets.use = c(8,10),  title.name = paste0("Signaling received by Inflam.DC and .TC - ", names(object.list)[i]), legend.pos.x = 10)
# }

##可视化，表达情况
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("left", "right")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T, type = "violin")


####
## 美化
p <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1,2,3,4,6,7,8), remove.isolate = FALSE,pairLR.use = pairLR.use, grid.on=T,color.grid = "black")

# 看颜色范围
range(p$data$prob,na.rm = T)
summary(p$data$prob,na.rm = T)

p1 <- p + scale_size_continuous(range = c(4, 8), guide = "none") +   # 调整气泡大小范围
  scale_color_gradientn(
    colours = c("#2760a9", "white", "#e50f20"),  # 定义颜色向量
    values = scales::rescale(c(0.2, 0.35, 0.5)),  # 定义颜色映射的范围
    name = "Commun. Prob." ) +  # 图例标题 
  xlab(label = NULL) + 
  ylab(label = NULL) + 
  geom_vline(xintercept = seq(1.5, length(unique(df.net$source)) - 0.5, 1)[1:6],lwd = 0.5) + ## 根据 netVisual_bubble 函数的源码，修改格子线的粗细
  geom_hline(yintercept = seq(1.5, length(unique(df.net$interaction_name_2)) - 0.5, 1)[1:3], lwd = 0.5) + 
  theme(axis.title.x = element_text(size = 16),  # 设置 x 轴标题字体大小
        axis.title.y = element_text(size = 16),  # 设置 y 轴标题字体大小
        axis.text.x = element_text(size = 14),  # 设置 x 轴刻度标签字体大小
        axis.text.y = element_text(size = 14, face = "italic"),   # 设置 y 轴刻度标签字体大小
        panel.border = element_rect(color = "black", fill=NA, size=1),  # 设置四周边框的颜色和粗细
        legend.key.size = unit(0.6, 'cm'),  # 设置图例键的大小
        legend.text = element_text(size = 12),  # 设置图例文本的大小
        legend.title = element_text(size = 13)
  ) 
p1

# 保存
ggsave(filename = "cellchat_bubble.pdf", width = 7.6, height = 4.5,plot = p1)













