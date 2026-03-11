#############after scVI
import multiprocessing
num_cores = multiprocessing.cpu_count()
workers = max(1, num_cores - 10)  # 留一个核心给系统 -10
print(f"Using {workers} workers for multiprocessing")


import omicverse as ov
print(f"omiverse version:{ov.__version__}")
import scanpy as sc
print(f"scanpy version:{sc.__version__}")
import scvi
print(f"scvi version:{scvi.__version__}")
import scib
print(f"scib version:{scib.__version__}")
ov.ov_plot_set()
import os
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42  # enables correct plotting of text for PDFs
os.chdir('/public8/lilab/student/htang/Metastasis/data/DISCO/0_rawdata/1_anno/norm/finished/decontx/merge/results')
adata=sc.read_h5ad('./pancancer_scvi_integrated_new112_new.h5ad')
adata = adata[:, adata.var.highly_variable]
adata
adata.X.max()

sc.tl.leiden(adata, key_added="leiden_res0.01", resolution=0.01, neighbors_key="scVI")
sc.tl.leiden(adata, key_added="leiden_res0.03", resolution=0.03, neighbors_key="scVI")
sc.tl.leiden(adata, key_added="leiden_res0.05", resolution=0.05, neighbors_key="scVI")
sc.tl.leiden(adata, key_added="leiden_res0.1", resolution=0.1, neighbors_key="scVI")
sc.tl.leiden(adata, key_added="leiden_res0.2", resolution=0.2, neighbors_key="scVI")
sc.tl.leiden(adata, key_added="leiden_res0.3", resolution=0.3, neighbors_key="scVI")
sc.tl.leiden(adata, key_added="leiden_res0.4", resolution=0.4, neighbors_key="scVI")
sc.tl.leiden(adata, key_added="leiden_res0.5", resolution=0.5, neighbors_key="scVI")
sc.tl.leiden(adata, key_added="leiden_res0.8", resolution=0.8, neighbors_key="scVI")
sc.tl.leiden(adata, key_added="leiden_res1.0", resolution=1.0, neighbors_key="scVI")

#######并行
from joblib import Parallel, delayed
import scanpy as sc
import anndata as ad
# 拷贝原始数据用于并行安全处理
import copy
adata_base = adata.copy()
# 并行函数
def run_leiden(res):
    adata_copy = adata_base.copy()
    key = f"leiden_res{res}"
    sc.tl.leiden(adata_copy, key_added=key, resolution=res, neighbors_key="scVI")
    return key, adata_copy.obs[key]

# 分辨率列表
res_list = [0.05, 0.1, 0.2, 0.3, 0.4]

# 并行执行
results = Parallel(n_jobs=5)(delayed(run_leiden)(r) for r in res_list)

# 将结果写回原始 adata（主线程安全操作）
for key, leiden_series in results:
    adata.obs[key] = leiden_series
#######################

sc.pl.umap(adata, neighbors_key='scVI',color=["leiden_res0.2"])
sc.pl.umap(adata, color=['leiden_res0.2'], legend_loc='on data', legend_fontsize=10, title="", frameon=False,
                  )
sc.pl.umap(adata, color=['leiden_res0.4'], legend_loc='on data', legend_fontsize=10, title="", frameon=False,
                   save='celltype_leiden_res0.4.png')
adata.write_h5ad('/public8/lilab/student/htang/Metastasis/data/DISCO/0_rawdata/1_anno/norm/finished/decontx/anno/pancancer_res.h5ad')
sc.tl.rank_genes_groups(adata, 'leiden_res0.2', method='wilcoxon')
import matplotlib.pyplot as plt
plt.rcParams['figure.subplot.bottom'] = 0.2  # 设置底部边距
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, fontsize=10)    ##foldchange=1.25

# 绘图
ax = sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, fontsize=10, show=False)
# 保存图像
plt.savefig("rank_genes_groups_top20.png", dpi=300, bbox_inches='tight')
# 如果需要，显示图形（可选）
# plt.show()

adata.obs['raw_celltype'].to_csv('brain_metastasis_raw_celltype.csv', header=True)
adata.write_h5ad('pancancer_res0.2_DEG.h5ad')
adata=sc.read_h5ad('pancancer_res0.2_DEG.h5ad')

#提取每组前50个基因
import pandas as pd
import numpy as np

# 将差异分析结果转换为 pandas 格式
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names  # 所有 cluster 名

# 构建 DataFrame
ranked_genes = {}
for group in groups:
    # 每个group的gene名、score
    names = result['names'][group]
    scores = result['scores'][group]

    # 组合为 DataFrame，按 score 降序排序
    df = pd.DataFrame({
        'gene': names,
        'score': scores
    }).sort_values('score', ascending=False).head(50)  # 取前20

    ranked_genes[group] = df

# 打印所有group前20基因
for group, df in ranked_genes.items():
    print(f"\nGroup {group} Top 50 genes by score:")
    print(df.to_string(index=False))

import pandas as pd

# 假设 ranked_genes 是一个 dict：{cluster: DataFrame}
# 每个 DataFrame 包含 columns 如 ['gene', 'score', 'logfoldchanges', 'pvals', 'pvals_adj']

all_dfs = []

for group, df in ranked_genes.items():
    df_top50 = df.head(50).copy()  # 前20个基因
    df_top50["cluster"] = group  # 添加cluster列
    all_dfs.append(df_top50)

# 合并所有数据
combined_df = pd.concat(all_dfs, ignore_index=True)

# 调整列顺序
cols = ['cluster'] + [col for col in combined_df.columns if col != 'cluster']
combined_df = combined_df[cols]

# 保存为 CSV 文件
combined_df.to_csv("ranked_genes_top50_by_score.csv", index=False)


###marker
genes = {
        'Epithelial cells': ["EPCAM",'KRT19','KRT15','KRT17','FABP1','FBP1','FOXJ1','WFDC2','PARD3'],
        'T and NK cells': ["PTPRC",'CD3D', 'CD3E', 'CD3G','CD2','CD4','CD8A','NKG7', 'GNLY'],
        'B cells': ["CD79A", "CD79B", "MS4A1", "CD19"],
        'Plasma cells': ['CD79A', 'JCHAIN', 'MZB1', 'IGHG1'],
        'Mast cells':['CST3', 'KIT', 'TPSB2', 'TPSAB1', 'MS4A2'],
        'Dendritic cells':['LILRA4','CXCR3','IRF7','IL3RA'],
        'Neutrophils':['CXCR2','CSF3R','FCGR3B'],
        'Endothelial cells': ['PECAM1', 'CLDN5', 'VWF'],
        'Fibroblasts': ["COL1A1", "DCN", "COL1A2", "LUM", "C1S"],  # ACTA2; THY1
        'Monocytes and Macrophages': ['CD68', 'CD163', 'CD14',"VCAN"],
        'Neurons and Glial cells': ['SYT1','NCAM1','APLP1','PLP1','S100B'],
        'Melanocytes': ['MLANA','MITF','PMEL','LRMDA'],
}
sc.pl.dotplot(adata, genes, groupby='leiden_res0.2', swap_axes=True)

import matplotlib.pyplot as plt
sc.pl.dotplot(adata, genes, groupby='leiden_res0.2', swap_axes=True, save='dotplot_anno01_cell.png')
plt.show()


# 细胞注释
# 定义映射，将多个簇合并为一个
cluster_mapping = {
    '0': 'NK and T cells',
    '1': 'Epithelial cells',
    '2': 'Monocytes and Macrophages',   
    '3': 'Fibroblasts',
    '4': 'Fibroblasts',
    '5': 'Melanocytes',
    '6': 'B cells',
    '7': 'Epithelial cells',
    '8': 'Melanocytes',
    '9': 'NK and T cells',
    '10': 'Endothelial cells',
    '11': 'Plasma cells',
    '12': 'Neutrophils',     
    '13': 'Epithelial cells',   
    '14': 'Melanocytes',    
    '15': 'Neurons and Glial cells',      
    '16': 'Mast cells',
    '17': 'Monocytes and Macrophages', 
    '18': 'Dendritic cells',
    '19': 'Epithelial cells',   
    '20': 'Epithelial cells',     
    '21': 'Epithelial cells',    
    '22': 'NK and T cells',    
    '23': 'NK and T cells',  
    '24': 'Epithelial cells',
    '25': 'Monocytes and Macrophages',
    '26': 'Endothelial cells',
    '27': 'NK and T cells',    
    '28': 'Plasma cells',
    '29': 'Epithelial cells',
    '30': 'Monocytes and Macrophages'    
    # 继续映射其他簇
}

# 应用映射
adata.obs['anno01'] = adata.obs['leiden_res0.2'].map(cluster_mapping)
# 注释后绘图
sc.pl.umap(
    adata, color="anno01", legend_loc="on data", title="", frameon=False, save=".pdf"
)
sc.pl.umap(
    adata, color="anno01", legend_loc="on data", title="", frameon=False, save=".png"
)
sc.pl.umap(
    adata, color="anno01", legend_loc="right margin", title="", frameon=False, save="2.pdf"
)

adata.write_h5ad('pancancer_res0.2_anno01.h5ad')

##########NKT
adata_NKT  = raw_adata[raw_adata.obs['anno01']=='NK and T cells'].copy()
adata_NKT.write_h5ad('/public8/lilab/student/htang/Metastasis/data/DISCO/0_rawdata/1_anno/norm/finished/decontx/anno/pancancer_res0.2_anno01_NKT.h5ad')
import omicverse as ov
print(f"omiverse version:{ov.__version__}")
import scanpy as sc
print(f"scanpy version:{sc.__version__}")
import scvi
print(f"scvi version:{scvi.__version__}")

ov.ov_plot_set()
import os
os.chdir('/public8/lilab/student/htang/Metastasis/data/DISCO/0_rawdata/1_anno/norm/finished/decontx/anno')
########NK and T cells
adata_NKT=sc.read_h5ad('./pancancer_res0.2_anno01_counts_NKT.h5ad')
adata_NKT.layers['counts'].max()
sc.pp.normalize_total(adata_NKT, target_sum=1e4)
sc.pp.log1p(adata_NKT)
adata_NKT.raw = adata_NKT
sc.pp.highly_variable_genes(adata_NKT,n_top_genes=3000, flavor='seurat')
adata_NKT.X.max()
sc.pp.scale(adata_NKT, max_value=10)
sc.pp.pca(adata_NKT)
#sc.tl.pca(adata_NKT, n_comps = 30,use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata_NKT,n_pcs = 15)
#sc.tl.leiden(adata_NKT, resolution=0.3, key_added='leiden_res0.3')
sc.tl.leiden(adata_NKT,flavor="igraph",n_iterations=2,resolution=0.5,key_added='leiden_NKT_res0.5')
sc.tl.umap(adata_NKT)
sc.pl.umap(adata_NKT,color='leiden_NKT_res0.5')

sc.pl.umap(adata_NKT, color=['leiden_NKT_res0.5'], legend_loc='on data',legend_fontsize=10, title="", frameon=False)
sc.pl.umap(adata_NKT, color=['leiden_NKT_res0.5'], legend_fontsize=10, title="", frameon=False,save='pancancer_anno02_NKT_celltype_0.5_1.png')
sc.tl.rank_genes_groups(adata_NKT, 'leiden_NKT_res0.5', method='wilcoxon')
#sc.tl.rank_genes_groups(adata_End, 'leiden1', method='wilcoxon')
import matplotlib.pyplot as plt
plt.rcParams['figure.subplot.bottom'] = 0.2  # 设置底部边距
sc.pl.rank_genes_groups(adata_NKT, n_genes=20, sharey=False, fontsize=10)
adata_NKT.write("./pancancer_res0.2_anno01_counts_NKT_res0.5_DEG.h5ad")
adata_NKT=sc.read_h5ad("./pancancer_res0.2_anno01_counts_NKT_res0.5_DEG.h5ad")
#dotplot
sc.pl.dotplot(adata, marker_genes, groupby="leiden0.5");

#绘制小提琴图
sc.pl.stacked_violin(adata, marker_genes, groupby="leiden");


# 绘图
ax = sc.pl.rank_genes_groups(adata_NKT, n_genes=20, sharey=False, fontsize=10, show=False)
# 保存图像
plt.savefig("NKT_rank_genes_groups_top20.png", dpi=300, bbox_inches='tight')
# 如果需要，显示图形（可选）
# plt.show()


#提取每组前50个基因
import pandas as pd
import numpy as np

# 将差异分析结果转换为 pandas 格式
result = adata_NKT.uns['rank_genes_groups']
groups = result['names'].dtype.names  # 所有 cluster 名

# 构建 DataFrame
ranked_genes = {}
for group in groups:
    # 每个group的gene名、score
    names = result['names'][group]
    scores = result['scores'][group]

    # 组合为 DataFrame，按 score 降序排序
    df = pd.DataFrame({
        'gene': names,
        'score': scores
    }).sort_values('score', ascending=False).head(50)  # 取前20

    ranked_genes[group] = df

# 打印所有group前20基因
for group, df in ranked_genes.items():
    print(f"\nGroup {group} Top 50 genes by score:")
    print(df.to_string(index=False))

import pandas as pd

# 假设 ranked_genes 是一个 dict：{cluster: DataFrame}
# 每个 DataFrame 包含 columns 如 ['gene', 'score', 'logfoldchanges', 'pvals', 'pvals_adj']

all_dfs = []

for group, df in ranked_genes.items():
    df_top50 = df.head(50).copy()  # 前20个基因
    df_top50["cluster"] = group  # 添加cluster列
    all_dfs.append(df_top50)

# 合并所有数据
combined_df = pd.concat(all_dfs, ignore_index=True)

# 调整列顺序
cols = ['cluster'] + [col for col in combined_df.columns if col != 'cluster']
combined_df = combined_df[cols]

# 保存为 CSV 文件
combined_df.to_csv("NKT_ranked_genes_top50_by_score.csv", index=False)





#########################
genes_NKTcell = {
        'CD4+ cell': ["CD3D", "CD3E", "CD3G",'CD4'],
        'CD8+ cell': ['GZMK','CD8A','CD8B'],
        'NK cell': ["KLRD1", "NKG7", "GNLY"],
    }

import matplotlib.pyplot as plt
sc.pl.dotplot(adata_NKT, genes_NKTcell, groupby='leiden_NKT_res0.5', swap_axes=True, save='dotplot_NKTcell.png')
plt.show()


# 细胞注释
# 定义映射，将多个簇合并为一个
cluster_mapping = {
    '0': 'CD8+ T cells',
    '1': 'NK cells',
    '2': 'CD8+ T cells',   ##???CD8+ T cells
    '3': 'CD4+ T cells',   #Treg
    '4': 'CD4+ T cells',
    '5': 'CD8+ T cells',
    '6': 'CD8+ T cells',
    '7': 'CD8+ T cells',
    '8': 'NK cells',
    '9': 'CD8+ T cells',
    '10': 'CD8+ T cells',
    '11': 'NK cells',
    '12': 'CD4+ T cells',    ##Monocytes and macrophages?
    '13': 'CD8+ T cells',   #???
    '14': 'CD4+ T cells',    ###???
    '15': 'CD4+ T cells',     #Mature Neurons, Oligodendrocytes, Astrocytes
    '16': 'CD8+ T cells',
    '17': 'CD8+ T cells',#Osteoclasts, Macrophages, Microglia
    '18': 'CD8+ T cells',
    '19': 'CD8+ T cells',  #alveolar epithelial cells
    '20': 'CD4+ T cells',    #multiciliated epithelial cells
    '21': 'CD8+ T cells',    #Hepatocytes
    '22': 'CD8+ T cells',    #Hepatocytes
}

# 应用映射
adata_NKT.obs['anno02_NKT'] = adata_NKT.obs['leiden_NKT_res0.5'].map(cluster_mapping)
adata_NKT.obs.rename(columns={'anno02_ NKT': 'anno02_NKT'}, inplace=True)

# 注释后绘图
sc.pl.umap(
    adata_NKT, color="anno02_NKT", legend_loc="on data", title="", frameon=False, save=".anno02_ NKT.pdf"
)
sc.pl.umap(
    adata_NKT, color="anno02_NKT", legend_loc="on data", title="", frameon=False, save=".anno02_ NKT.png"
)
sc.pl.umap(
    adata_NKT, color="anno02_NKT", legend_loc="right margin", title="", frameon=False, save=".anno02_ NKT2.png"
)

adata_NKT.write_h5ad('pancancer_res0.2_anno01_NKT_anno02_NKT.h5ad')


sc.pl.dotplot(adata_NKT, genes_NKTcell,'anno02_NKT', dendrogram=True,save='anno02_NKT_5.png')

###################MonMac
adata_MonMac = raw_adata[raw_adata.obs['anno01']=='Monocytes and Macrophages'].copy()
adata_MonMac.write_h5ad('/public8/lilab/student/htang/Metastasis/data/DISCO/0_rawdata/1_anno/norm/finished/decontx/anno/pancancer_res0.2_anno01_Mon_Mac.h5ad')
adata_MonMac=sc.read_h5ad('./pancancer_res0.2_anno01_Mon_Mac.h5ad')
adata_MonMac.layers['counts'].max()
sc.pp.normalize_total(adata_MonMac, target_sum=1e4)
sc.pp.log1p(adata_MonMac)
adata_MonMac.raw = adata_MonMac
sc.pp.highly_variable_genes(adata_MonMac,n_top_genes=3000, flavor='seurat')
adata_MonMac.X.max()
sc.pp.scale(adata_MonMac, max_value=10)
sc.pp.pca(adata_MonMac)
#sc.tl.pca(adata_MonMac, n_comps = 30,use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata_MonMac,n_pcs = 15)
#sc.tl.leiden(adata_MonMac, resolution=0.3, key_added='leiden_res0.3')
sc.tl.leiden(adata_MonMac,flavor="igraph",n_iterations=2,resolution=0.5,key_added='leiden_MonMac_res0.5')
#sc.tl.umap(adata_MonMac)
sc.pl.umap(adata_MonMac,color='leiden_MonMac_res0.5')

sc.pl.umap(adata_MonMac, color=['leiden_MonMac_res0.5'], legend_loc='on data',legend_fontsize=10, title="", frameon=False)
sc.pl.umap(adata_MonMac, color=['leiden_MonMac_res0.5'], legend_fontsize=10, title="", frameon=False,save='pancancer_anno02_MonMac_celltype_0.5_1.png')
sc.tl.rank_genes_groups(adata_MonMac, 'leiden_MonMac_res0.5', method='wilcoxon')
#sc.tl.rank_genes_groups(adata_End, 'leiden1', method='wilcoxon')
import matplotlib.pyplot as plt
plt.rcParams['figure.subplot.bottom'] = 0.2  # 设置底部边距
sc.pl.rank_genes_groups(adata_MonMac, n_genes=20, sharey=False, fontsize=10)
adata_MonMac.write("./pancancer_res0.2_anno01_counts_MonMac_res0.5_DEG.h5ad")
adata_MonMac=sc.read_h5ad("./pancancer_res0.2_anno01_counts_MonMac_res0.5_DEG.h5ad")

# 绘图
ax = sc.pl.rank_genes_groups(adata_MonMac, n_genes=20, sharey=False, fontsize=10, show=False)
# 保存图像
plt.savefig("NKT_rank_genes_groups_top20.png", dpi=300, bbox_inches='tight')
# 如果需要，显示图形（可选）
# plt.show()


#提取每组前50个基因
import pandas as pd
import numpy as np

# 将差异分析结果转换为 pandas 格式
result = adata_MonMac.uns['rank_genes_groups']
groups = result['names'].dtype.names  # 所有 cluster 名

# 构建 DataFrame
ranked_genes = {}
for group in groups:
    # 每个group的gene名、score
    names = result['names'][group]
    scores = result['scores'][group]

    # 组合为 DataFrame，按 score 降序排序
    df = pd.DataFrame({
        'gene': names,
        'score': scores
    }).sort_values('score', ascending=False).head(50)  # 取前20

    ranked_genes[group] = df

# 打印所有group前20基因
for group, df in ranked_genes.items():
    print(f"\nGroup {group} Top 50 genes by score:")
    print(df.to_string(index=False))

import pandas as pd

# 假设 ranked_genes 是一个 dict：{cluster: DataFrame}
# 每个 DataFrame 包含 columns 如 ['gene', 'score', 'logfoldchanges', 'pvals', 'pvals_adj']

all_dfs = []

for group, df in ranked_genes.items():
    df_top50 = df.head(50).copy()  # 前20个基因
    df_top50["cluster"] = group  # 添加cluster列
    all_dfs.append(df_top50)

# 合并所有数据
combined_df = pd.concat(all_dfs, ignore_index=True)

# 调整列顺序
cols = ['cluster'] + [col for col in combined_df.columns if col != 'cluster']
combined_df = combined_df[cols]

# 保存为 CSV 文件
combined_df.to_csv("MonMac_ranked_genes_top50_by_score.csv", index=False)





#########################
genes_MonMac = {
        'Macrophages': ["CD163", "CD68",'C1QA','C1QB','C1QC','MMP9'],
        'Monocytes': ['FCN1','VCAN','APOBEC3A','THBS1']
    }


import matplotlib.pyplot as plt
sc.pl.dotplot(adata_MonMac, genes_MonMac, groupby='leiden_MonMac_res0.5', swap_axes=True, save='dotplot_MonMac_res0.5.png')
plt.show()


# 细胞注释
# 定义映射，将多个簇合并为一个
cluster_mapping = {
    '0': 'Macrophages',
    '1': 'Macrophages',
    '2': 'Macrophages',   ##???CD8+ T cells
    '3': 'Macrophages',   #Treg
    '4': 'Macrophages',
    '5': 'Monocytes',
    '6': 'Macrophages',
    '7': 'Macrophages',
    '8': 'Macrophages',
    '9': 'Macrophages',
    '10': 'Macrophages',
    '11': 'Macrophages',
    '12': 'Macrophages',    ##Monocytes and macrophages?
    '13': 'Macrophages',   #???
    '14': 'Macrophages',    ###???
    '15': 'Monocytes',     #Mature Neurons, Oligodendrocytes, Astrocytes
    '16': 'Macrophages',
    '17': 'Monocytes',#Osteoclasts, Macrophages, Microglia
    '18': 'Monocytes',
    '19': 'Macrophages',  #alveolar epithelial cells
    '20': 'Macrophages',    #multiciliated epithelial cells
    '21': 'Macrophages',    #Hepatocytes
    '22': 'Macrophages',    #Hepatocytes
    '23': 'Macrophages',  # Hepatocytes
    '24': 'Macrophages',  # Hepatocytes
}

# 应用映射
adata_MonMac.obs['anno02_MonMac'] = adata_MonMac.obs['leiden_MonMac_res0.5'].map(cluster_mapping)

# 注释后绘图
sc.pl.umap(
    adata_MonMac, color="anno02_MonMac", legend_loc="on data", title="", frameon=False, save=".anno02_MonMac.pdf"
)
sc.pl.umap(
    adata_MonMac, color="anno02_MonMac", legend_loc="on data", title="", frameon=False, save=".anno02_MonMac.png"
)
sc.pl.umap(
    adata_MonMac, color="anno02_MonMac", legend_loc="right margin", title="", frameon=False, save=".anno02_MonMac_2.png"
)

adata_MonMac.write_h5ad('pancancer_res0.2_anno01_MonMac_anno02_MonMac.h5ad')
