# ---****----utf-8---****----
# @File  : tmp.py
# @Author: Tang Hai
# @email : tangh25@mail2.sysu.edu.cn
# @Date  :  2025/10/23
import os
import pandas as pd
import numpy as np
import scanpy as sc
from scib_metrics.benchmark import Benchmarker
from scib_metrics.benchmark import Benchmarker, BioConservation
from scib_metrics.benchmark import BatchCorrection
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42  # enables correct plotting of text for PDFs
print('20250208_13w')
os.chdir('/public8/lilab/student/htang/Metastasis/analysis/data/dataprocess/merge/results')
adata=sc.read_h5ad('./panMET_pca_harmony_BBKNN_combat_scanorama_scVI_anno_celltypist.h5ad')
type(adata.X)
# 2. 查看结果
print(f"共筛选到 {adata.var['highly_variable'].sum()} 个高变基因")
# 3. 提取高变基因矩阵
#adata = adata[:, adata.var['highly_variable']].copy()
SEED = 42
np.random.seed(SEED)
sc.pp.subsample(adata, n_obs=100_000)
# 抽取 1/20 (即 5%) 的细胞
#sc.pp.subsample(adata, fraction=0.05)
adata.var["highly_variable"] = np.asarray(adata.var["highly_variable"].astype(bool))
# import scipy.sparse as sp
# adata.X = sp.csr_matrix(adata.X)
# print("✅ 已将 adata.X 转换为稀疏矩阵 (CSR format)")
adata.obs['majority_voting'].unique()
sc.tl.pca(adata)
#######################################################

#使用scib包来对不同去批次效应方法进行测评
from scib_metrics.benchmark import Benchmarker
from scib_metrics.benchmark import Benchmarker, BioConservation
from scib_metrics.benchmark import BatchCorrection
# biocons = BioConservation(isolated_labels=False, silhouette_label= False)
# batchcor = BatchCorrection(kbet_per_label=False)

bm = Benchmarker(
    adata,
    batch_key="sample",
    label_key="majority_voting",
    embedding_obsm_keys=["X_pca", "X_harmony", "X_scanorama", "X_scVI", "X_BBKNN"],   #, "X_combat"
    n_jobs=16,    #8
)
bm.benchmark()
df=bm.get_results(min_max_scale=False)
print(df)
df.to_csv("batch_integration_benchmark_results0208_13w.csv", index=True)
print("✅ Benchmark 结果已保存为 batch_integration_benchmark_results0208_13w.csv")
#
# bm.plot_results_table(min_max_scale=False)
# fig = bm.plot_results_table(min_max_scale=False)
# 获取表格背后的 matplotlib figure 对象
import matplotlib.pyplot as plt
# 绘制表格
fig_table = bm.plot_results_table(min_max_scale=False)  # 返回的是 Table 对象
# 获取当前 Figure
fig = plt.gcf()
fig.savefig("benchmark_results_table0208_13w.png", dpi=300, bbox_inches='tight')
fig.savefig("benchmark_results_table0208_13w.pdf", dpi=300, bbox_inches='tight')

plt.close(fig)


adata.write_h5ad('./panMET_pca_harmony_BBKNN_combat_scanorama_scVI_anno_celltypist_benchmark_0208_13w.h5ad')
print("✅ Benchmark 结果图已保存为 benchmark_results_table0208_13w.png")


