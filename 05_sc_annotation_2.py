############## sub cell type annotation using TOSICA #############################################

# Import packages
from TOSICA import train
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "1"

import TOSICA
import scanpy as sc
import numpy as np
import warnings
warnings.filterwarnings ("ignore")

import torch
print(torch.__version__)
print(torch.cuda.get_device_capability(device=None), torch.cuda.get_device_name(device=None))

# Set working directory
os.chdir('./TOSICA-main/')


############################
# Reference dataset preprocessing
############################

# Load reference fibroblast dataset
ref_adata = sc.read('./ref_adata.h5ad')

print(f"The max of ref_adata is {ref_adata.X.max()}")

# Extract cell type labels
ref_adata.obs["ref_celltype"] = ref_adata.obs["MajorCluster"].str.split("_", n=1).str[1]
print(ref_adata.obs["ref_celltype"].unique())


############################
# Subsample reference cells (>30k)
############################

import random
random.seed(42)

cell_idx = list(random.sample(ref_adata.obs.index.tolist(), 30000))
ref_adata = ref_adata[cell_idx]

print(ref_adata.obs["ref_celltype"].value_counts())


############################
# Reverse log1p transformation
############################

ref_adata.X = np.expm1(ref_adata.X)
print(ref_adata.X.max())

# Normalize and log transform again
ref_adata.var_names_make_unique()

ref_adata.layers["counts"] = ref_adata.X.copy()

sc.pp.normalize_total(ref_adata, target_sum=1e4)
sc.pp.log1p(ref_adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(
    ref_adata,
    n_top_genes=5000,
)


############################
# Query dataset preprocessing
############################

# Load fibroblast cells from pancancer dataset
adata=sc.read_h5ad('./counts_major_cell_type.h5ad')

adata.layers['counts']=adata.X.copy()

adata.var_names_make_unique()

adata.layers["counts"] = adata.X.copy()

# Normalize and log transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata)


############################
# Identify shared genes
############################

ref_adata.var_names_make_unique()
adata.var_names_make_unique()

# Keep HVGs in reference
ref_adata = ref_adata[:, ref_adata.var.highly_variable]

# Intersect genes between reference and query
ret_gene=list(set(adata.var_names) & set(ref_adata.var_names))
print(len(ret_gene))

# Subset both datasets to shared genes
adata=adata[:,ret_gene]
ref_adata=ref_adata[:,ret_gene]

# Save processed datasets
ref_adata.write_h5ad('.ref_celltype_50000unscale.h5ad')
adata.write_h5ad('./counts_major_cell_type_5000unscale.h5ad')

print(f"The max of ref_adata is {ref_adata.X.max()}, query_data is {adata.X.max()}")

############################
# TOSICA model training
############################

import time

ref_adata = ref_adata[:,ref_adata.var_names]

print(ref_adata)
print(ref_adata.obs.ref_celltype.value_counts())

t_start = time.time()

TOSICA.train(
    ref_adata,
    gmt_path='human_gobp',
    label_name='ref_celltype',
    batch_size=256,
    epochs=50,
    depth=1,
    project='hGOBP_demo_Mac_epo50_5'
)

t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes!")


############################
# Prediction using trained model
############################
model_weight_path = './model.pth'

new_adata = TOSICA.pre(
    adata,
    model_weight_path = model_weight_path,
    project='epo100_depth',
    depth=2
)

print(new_adata.obs['Prediction'].unique().tolist())
print(new_adata.obs['Prediction'].value_counts())


############################
# Filter low-frequency predicted cell types
############################

cell_idx=new_adata.obs['Prediction'].value_counts()[new_adata.obs['Prediction'].value_counts()<200].index
new_adata=new_adata[~new_adata.obs['Prediction'].isin(cell_idx)]

print(new_adata.obs['Prediction'].value_counts())


############################
# Evaluate prediction confidence
############################

pred_values = new_adata.obs['Probability'].astype(float)

low_prediction_count = (pred_values < 0.3).sum()
print("Number of cells with probability < 0.3:", low_prediction_count)

low_pred_cells = new_adata.obs[pred_values < 0.5]
low_pred_cells.to_csv("low_confidence_cells.csv")


############################
# Identify pathway attention signals
############################

sc.tl.rank_genes_groups(new_adata, 'Prediction', method='wilcoxon')

sc.pl.rank_genes_groups_dotplot(
    new_adata,
    n_genes=3,
    standard_scale='var',
    save='cell_gene_subanno_epo.png'
)


############################
# Convert metadata types
############################

for col in new_adata.obs.columns:
    if new_adata.obs[col].dtype == 'object':
        try:
            new_adata.obs[col] = pd.to_numeric(new_adata.obs[col], errors='raise')
        except:
            new_adata.obs[col] = new_adata.obs[col].astype(str)


############################
# Save annotated dataset
############################

new_adata.write_h5ad('./TOSICA_anno/Pcancer_anno01_anno2_epo_TOSICA_anno.h5ad')

# Reload original dataset if needed
adata=sc.read_h5ad('./pancancer_anno02_counts.h5ad')
