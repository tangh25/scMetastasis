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
        'Epithelial cells': ["EPCAM",'KRT19','KRT15','KRT17'],
        'T and NK cells': ["PTPRC",'CD3D', 'CD3E', 'CD3G','CD2','CD4','CD8A','NKG7', 'GNLY'],
        'B cells': ["CD79A", "CD79B", "MS4A1", "CD19"],
        'Plasma cells': ['CD79A', 'JCHAIN', 'MZB1', 'IGHG1'],
        'Mast cells':['CST3', 'KIT', 'TPSB2', 'TPSAB1', 'MS4A2'],
        'Dendritic cells':['LILRA4','CXCR3','IRF7','IL3RA'],
        'Neutrophils':['CXCR2','CSF3R','FCGR3B'],
        'Endothelial cells': ['PECAM1', 'CLDN5', 'VWF'],
        'Fibroblasts': ["COL1A1", "DCN", "COL1A2", "LUM", "C1S"], 
        'Monocytes and Macrophages': ['CD68', 'CD163', 'CD14',"VCAN"],
        'Neurons and Glial cells': ['SYT1','NCAM1','APLP1','PLP1','S100B'],
        'Melanoma cells': ['MLANA','MITF','PMEL','LRMDA'],
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
    '2': 'Monocytes and Macrophages',   #Macrophages
    '3': 'Fibroblasts',
    '4': 'Fibroblasts',
    '5': 'Melanocytes',
    '6': 'B cells',
    '7': 'Epithelial cells',
    '8': 'Melanocytes',
    '9': 'NK and T cells',
    '10': 'Endothelial cells',
    '11': 'Plasma cells',
    '12': 'Neutrophils',    ##Monocytes and macrophages?
    '13': 'Epithelial cells',   #???
    '14': 'Melanocytes',    ###???
    '15': 'Neurons and Glial cells',     #Mature Neurons, Oligodendrocytes, Astrocytes
    '16': 'Mast cells',
    '17': 'Monocytes and Macrophages',#Osteoclasts, Macrophages, Microglia
    '18': 'Dendritic cells',
    '19': 'Epithelial cells',  #alveolar epithelial cells
    '20': 'Epithelial cells',    #multiciliated epithelial cells
    '21': 'Epithelial cells',    #Hepatocytes
    '22': 'NK and T cells',   #?
    '23': 'NK and T cells',  #?
    '24': 'Epithelial cells',
    '25': 'Monocytes and Macrophages',
    '26': 'Endothelial cells',
    '27': 'NK and T cells',    ###?
    '28': 'Plasma cells',
    '29': 'Epithelial cells',
    '30': 'Monocytes and Macrophages'    #Monocytes and macrophages
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
