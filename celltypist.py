# ---****----utf-8---****----
# @File  : PanMET_CellTypist.py
# @Author: Tang Hai
# @email : tangh25@mail2.sysu.edu.cn
# @Date  :  2025/07/22

##如何将细胞类型标签从 adata_James 转移到 adata_Elmentaite，并评估和可视化预测结果
import scanpy as sc
import celltypist
import time
import numpy as np

import os
os.chdir('/public8/lilab/student/htang/Metastasis/data/DISCO/0_rawdata/1_anno/norm/finished/decontx/anno/anno02/Celltypist/')
del adata
#logP  ref
adata_Elmentaite = sc.read('/public8/lilab/student/htang/Metastasis/data/DISCO/0_rawdata/1_anno/norm/finished/decontx/anno/anno02/ref/GSE154763_MajorClusters_Mac.h5ad')
print(adata_Elmentaite.X.max())
#单细胞转录组数据中的 NaN 值是可以（也通常应当）填充为 0
import numpy as np
adata_James.X = np.nan_to_num(adata_James.X, nan=0.0)
# 逆 log1p 转换：exp(x) - 1
import numpy as np
adata_Elmentaite.X = np.expm1(adata_Elmentaite.X)
print(adata_Elmentaite.X.max())
sc.pp.normalize_total(adata_Elmentaite, target_sum = 1e4)
sc.pp.log1p(adata_Elmentaite)
print(adata_Elmentaite.X.max())

### adata_James
adata_James = sc.read('/public8/lilab/student/htang/Metastasis/data/DISCO/0_rawdata/1_anno/norm/finished/decontx/anno/pancancer_res0.2_anno01_counts_Mac.h5ad')
print(adata_James.X.max())
sc.pp.normalize_total(adata_James, target_sum = 1e4)
sc.pp.log1p(adata_James)
print(adata_James.X.max())

print(adata_James.X.max())
print(adata_James.obs['MajorCluster'].unique())

# 134 cell types in the first data.
print(adata_Elmentaite.obs['MajorCluster'].unique().tolist())
# adata_Elmentaite.obs['Integrated_05'].unique()
# # 25 cell types in the second data.
# adata_James.obs.cell_type.unique()

#Transfer cell type labels from the first dataset to the second dataset
adata_Elmentaite.shape

# Sample 500 cells from each cell type for `adata_Elmentaite`.
# All cells from a given cell type will be selected if the cell type size is < 500.
#sampled_cell_index = celltypist.samples.downsample_adata(adata_Elmentaite, mode = 'each', n_cells = 500, by = 'Integrated_05', return_index = True)
##all cell
sampled_cell_index = adata_Elmentaite.obs_names.tolist()

print(f"Number of downsampled cells for training: {len(sampled_cell_index)}")

#(Suggested) Feature selection for the first dataset
# Use `celltypist.train` to quickly train a rough CellTypist model.
# You can also set `mini_batch = True` to enable mini-batch training.
t_start = time.time()
model_fs = celltypist.train(adata_Elmentaite[sampled_cell_index], 'MajorCluster', n_jobs = 30, max_iter = 5, use_SGD = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")
#根据与给定细胞类型相关的绝对回归系数，从每种细胞类型中绘制出排名前100的重要基因。对于只有几种细胞类型的数据显示集，您可能需要将顶级基因数量从100增加到例如300，以获得足够的基因用于最终使用。
gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -300, axis = 1)[:, -300:]
gene_index = np.unique(gene_index)
print(f"Number of genes selected: {len(gene_index)}")

#Model training and label transfer
# Add `check_expression = False` to bypass expression check with only a subset of genes.
t_start = time.time()
model = celltypist.train(adata_Elmentaite[sampled_cell_index, gene_index], 'MajorCluster', check_expression = False, n_jobs = 10, max_iter = 100)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes")

# Save the model.
model.write('celltypist_anno_model/model_from_pancancer_Mac.pkl')

#Next, use celltypist.annotate to predict adata_James using this model.
# CellTypist prediction without over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(adata_James, model = 'celltypist_anno_model/model_from_pancancer_Mac.pkl')
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")

# CellTypist prediction with over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(adata_James, model = 'celltypist_anno_model/model_from_pancancer_Mac.pkl', majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")

predictions.predicted_labels.head()


# You can also change the value of `use_as_prediction` to `predicted_labels` to compare the raw prediction result with the pre-defined cell types.
celltypist.dotplot(predictions, use_as_reference = 'Prediction', use_as_prediction = 'majority_voting')
celltypist.dotplot(predictions, use_as_reference = 'ann', use_as_prediction = 'majority_voting')

# Get an `AnnData` with predicted labels embedded into the cell metadata columns.
adata = predictions.to_adata()
adata.obs.iloc[:, -4:]
##查看 conf_score 小于 0.7 的细胞
adata.obs[adata.obs['conf_score'] < 0.7]
adata.obs[adata.obs['conf_score'] < 0.7]
(adata.obs['conf_score'] < 0.7).mean()

sc.pl.umap(adata, color = ['majority_voting'], legend_loc = 'on data')
sc.pl.umap(adata, color = ['Prediction', 'majority_voting'], legend_loc = 'on data')
adata.write_h5ad('panMET_celltypist_anno_Mac.h5ad')

adata.X.max()
adata.obs.index

import pandas as pd
anno2 =pd.read_csv('/public8/lilab/student/htang/Metastasis/data/DISCO/0_rawdata/1_anno/norm/finished/decontx/anno/anno02/Mac/Pcancer_anno01_Mac_cell_anno2_TOSICA_epo15.csv',index_col=0)
anno22 =anno2.iloc[:,0:2]
anno22.index = anno22.index.astype(str)
anno22.index

adata.obs = adata.obs.merge(anno22, left_index=True, right_index=True, how='left')
adata.obs=adata.obs.drop(columns=['Prediction_x','Probability_x'])
adata.obs=adata.obs.rename(columns={'Prediction_y':'Prediction','Probability_y':'Probability'})
#如果 'majority_voting' 和 'Prediction' 的值相同，就保留这个值；否则，设为 'unanno'；新结果保存在 adata.obs['final_label'] 中
import numpy as np
adata.obs['final_label'] = np.where(
    adata.obs['majority_voting'] == adata.obs['Prediction'],
    adata.obs['Prediction'],
    'Unknown'
)
adata.obs[['majority_voting', 'Prediction', 'final_label']].head()
adata.obs['final_label'].value_counts()
adata.obs.to_csv('pancacer_Mac_finallabel.csv')