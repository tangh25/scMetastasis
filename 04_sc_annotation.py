############# After scVI integration
# Detect number of CPU cores for multiprocessing
import multiprocessing
num_cores = multiprocessing.cpu_count()
workers = max(1, num_cores - 10)  # leave some cores for system
print(f"Using {workers} workers for multiprocessing")

# Import required libraries
import omicverse as ov
print(f"omiverse version:{ov.__version__}")

import scanpy as sc
print(f"scanpy version:{sc.__version__}")

# Set plotting style
ov.ov_plot_set()

import os
from matplotlib import rcParams

# Ensure correct font embedding in PDF
rcParams['pdf.fonttype'] = 42  

# Set working directory
os.chdir('./data')

# Load scVI-integrated dataset
adata = sc.read_h5ad('./pancancer_scvi_integrated.h5ad')

# Use only highly variable genes
adata = adata[:, adata.var.highly_variable]

# Check data
adata
adata.X.max()


############################
# Leiden clustering at multiple resolutions
############################

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


############################
# Parallel Leiden clustering
############################

from joblib import Parallel, delayed
import scanpy as sc
import anndata as ad

# Create a safe copy of dataset for parallel processing
import copy
adata_base = adata.copy()

# Define parallel Leiden function
def run_leiden(res):
    adata_copy = adata_base.copy()
    key = f"leiden_res{res}"
    sc.tl.leiden(adata_copy, key_added=key, resolution=res, neighbors_key="scVI")
    return key, adata_copy.obs[key]

# Resolution list
res_list = [0.05, 0.1, 0.2, 0.3, 0.4]

# Run clustering in parallel
results = Parallel(n_jobs=5)(delayed(run_leiden)(r) for r in res_list)

# Write results back to main AnnData object
for key, leiden_series in results:
    adata.obs[key] = leiden_series


############################
# UMAP visualization
############################

sc.pl.umap(adata, neighbors_key='scVI', color=["leiden_res0.2"])

sc.pl.umap(
    adata,
    color=['leiden_res0.2'],
    legend_loc='on data',
    legend_fontsize=10,
    title="",
    frameon=False
)

sc.pl.umap(
    adata,
    color=['leiden_res0.4'],
    legend_loc='on data',
    legend_fontsize=10,
    title="",
    frameon=False,
    save='celltype_leiden_res0.4.png'
)

# Save clustering result
adata.write_h5ad('./pancancer_res.h5ad')


############################
# Differential expression analysis
############################

# Identify marker genes for clusters
sc.tl.rank_genes_groups(adata, 'leiden_res0.2', method='wilcoxon')

import matplotlib.pyplot as plt
plt.rcParams['figure.subplot.bottom'] = 0.2

# Plot top marker genes
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, fontsize=10)

# Save figure
ax = sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, fontsize=10, show=False)
plt.savefig("rank_genes_groups_top20.png", dpi=300, bbox_inches='tight')


############################
# Extract top marker genes per cluster
############################

import pandas as pd
import numpy as np

# Get DEG results
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

ranked_genes = {}

# Extract top 50 genes per cluster
for group in groups:

    names = result['names'][group]
    scores = result['scores'][group]

    df = pd.DataFrame({
        'gene': names,
        'score': scores
    }).sort_values('score', ascending=False).head(50)

    ranked_genes[group] = df

# Print marker genes
for group, df in ranked_genes.items():
    print(f"\nGroup {group} Top 50 genes by score:")
    print(df.to_string(index=False))


############################
# Save DEG results
############################

all_dfs = []

for group, df in ranked_genes.items():

    df_top50 = df.head(50).copy()
    df_top50["cluster"] = group
    all_dfs.append(df_top50)

combined_df = pd.concat(all_dfs, ignore_index=True)

cols = ['cluster'] + [col for col in combined_df.columns if col != 'cluster']
combined_df = combined_df[cols]

combined_df.to_csv("ranked_genes_top50_by_score.csv", index=False)


############################
# Marker gene visualization
############################

genes = {
    'Epithelial cells': ["EPCAM",'KRT19','KRT15','KRT17'],
    'T and NK cells': ["PTPRC",'CD3D','CD3E','CD3G','CD2','CD4','CD8A','NKG7','GNLY'],
    'B cells': ["CD79A","CD79B","MS4A1","CD19"],
    'Plasma cells': ['CD79A','JCHAIN','MZB1','IGHG1'],
    'Mast cells':['CST3','KIT','TPSB2','TPSAB1','MS4A2'],
    'Dendritic cells':['LILRA4','CXCR3','IRF7','IL3RA'],
    'Neutrophils':['CXCR2','CSF3R','FCGR3B'],
    'Endothelial cells':['PECAM1','CLDN5','VWF'],
    'Fibroblasts':["COL1A1","DCN","COL1A2","LUM","C1S"],
    'Monocytes and Macrophages':['CD68','CD163','CD14',"VCAN"],
    'Neurons and Glial cells':['SYT1','NCAM1','APLP1','PLP1','S100B'],
    'Melanocytes':['MLANA','MITF','PMEL','LRMDA'],
}

# Dotplot of marker genes
sc.pl.dotplot(adata, genes, groupby='leiden_res0.2', swap_axes=True)

sc.pl.dotplot(
    adata,
    genes,
    groupby='leiden_res0.2',
    swap_axes=True,
    save='dotplot_anno01_cell.png'
)


############################
# Cell type annotation
############################

# Map clusters to cell types
cluster_mapping = {
    '0':'NK and T cells',
    '1':'Epithelial cells',
    '2':'Monocytes and Macrophages',
    '3':'Fibroblasts',
    '4':'Fibroblasts',
    '5':'Melanocytes',
    '6':'B cells',
    '7':'Epithelial cells',
    '8':'Melanocytes',
    '9':'NK and T cells',
    '10':'Endothelial cells',
    '11':'Plasma cells',
    '12':'Neutrophils',
    '13':'Epithelial cells',
    '14':'Melanocytes',
    '15':'Neurons and Glial cells',
    '16':'Mast cells',
    '17':'Monocytes and Macrophages',
    '18':'Dendritic cells',
    '19':'Epithelial cells',
    '20':'Epithelial cells',
    '21':'Epithelial cells',
    '22':'NK and T cells',
    '23':'NK and T cells',
    '24':'Epithelial cells',
    '25':'Monocytes and Macrophages',
    '26':'Endothelial cells',
    '27':'NK and T cells',
    '28':'Plasma cells',
    '29':'Epithelial cells',
    '30':'Monocytes and Macrophages'
}

# Apply annotation
adata.obs['anno01'] = adata.obs['leiden_res0.2'].map(cluster_mapping)

# Plot annotated UMAP
sc.pl.umap(adata, color="anno01", legend_loc="on data", frameon=False, save=".pdf")
sc.pl.umap(adata, color="anno01", legend_loc="on data", frameon=False, save=".png")
sc.pl.umap(adata, color="anno01", legend_loc="right margin", frameon=False, save="2.pdf")

# Save annotated dataset
adata.write_h5ad('pancancer_res0.2_anno01.h5ad')


################
########## Extract NK and T cells
# Subset NK and T cells from the annotated dataset
adata_NKT  = raw_adata[raw_adata.obs['anno01']=='NK and T cells'].copy()

# Save subset
adata_NKT.write_h5ad('./pancancer_res0.2_anno01_NKT.h5ad')


# Import packages and check versions
import scanpy as sc
print(f"scanpy version:{sc.__version__}")

ov.ov_plot_set()

import os
os.chdir('./data')


######## NK and T cells analysis
# Load NK and T cell dataset
adata_NKT=sc.read_h5ad('./pancancer_res0.2_anno01_counts_NKT.h5ad')

# Check counts layer
adata_NKT.layers['counts'].max()

# Normalize and log-transform
sc.pp.normalize_total(adata_NKT, target_sum=1e4)
sc.pp.log1p(adata_NKT)

# Store raw data
adata_NKT.raw = adata_NKT

# Identify highly variable genes
sc.pp.highly_variable_genes(adata_NKT,n_top_genes=3000, flavor='seurat')

adata_NKT.X.max()

# Scale data
sc.pp.scale(adata_NKT, max_value=10)

# PCA dimensionality reduction
sc.pp.pca(adata_NKT)

# Construct neighbor graph
sc.pp.neighbors(adata_NKT,n_pcs = 15)

# Leiden clustering
sc.tl.leiden(adata_NKT,flavor="igraph",n_iterations=2,resolution=0.5,key_added='leiden_NKT_res0.5')

# UMAP visualization
sc.tl.umap(adata_NKT)

sc.pl.umap(adata_NKT,color='leiden_NKT_res0.5')

# Plot UMAP with cluster labels
sc.pl.umap(adata_NKT, color=['leiden_NKT_res0.5'], legend_loc='on data',legend_fontsize=10, title="", frameon=False)

# Save UMAP
sc.pl.umap(adata_NKT, color=['leiden_NKT_res0.5'], legend_fontsize=10, title="", frameon=False,save='pancancer_anno02_NKT_celltype_0.5_1.png')


# Differential gene expression analysis
sc.tl.rank_genes_groups(adata_NKT, 'leiden_NKT_res0.5', method='wilcoxon')

import matplotlib.pyplot as plt

# Adjust figure margin
plt.rcParams['figure.subplot.bottom'] = 0.2

# Plot top marker genes
sc.pl.rank_genes_groups(adata_NKT, n_genes=20, sharey=False, fontsize=10)

# Save DEG results
adata_NKT.write("./pancancer_res0.2_anno01_counts_NKT_res0.5_DEG.h5ad")

# Reload dataset
adata_NKT=sc.read_h5ad("./pancancer_res0.2_anno01_counts_NKT_res0.5_DEG.h5ad")


# Dotplot visualization
sc.pl.dotplot(adata, marker_genes, groupby="leiden0.5")

# Violin plot visualization
sc.pl.stacked_violin(adata, marker_genes, groupby="leiden")


# Save DEG figure
ax = sc.pl.rank_genes_groups(adata_NKT, n_genes=20, sharey=False, fontsize=10, show=False)
plt.savefig("NKT_rank_genes_groups_top20.png", dpi=300, bbox_inches='tight')


# Extract top 50 marker genes for each cluster
import pandas as pd
import numpy as np

# Convert DEG results to pandas format
result = adata_NKT.uns['rank_genes_groups']
groups = result['names'].dtype.names

ranked_genes = {}

for group in groups:

    # Extract gene names and scores
    names = result['names'][group]
    scores = result['scores'][group]

    # Create dataframe and select top 50 genes
    df = pd.DataFrame({
        'gene': names,
        'score': scores
    }).sort_values('score', ascending=False).head(50)

    ranked_genes[group] = df


# Print top genes
for group, df in ranked_genes.items():
    print(f"\nGroup {group} Top 50 genes by score:")
    print(df.to_string(index=False))


# Combine results into a single table
all_dfs = []

for group, df in ranked_genes.items():

    df_top50 = df.head(50).copy()
    df_top50["cluster"] = group
    all_dfs.append(df_top50)

combined_df = pd.concat(all_dfs, ignore_index=True)

cols = ['cluster'] + [col for col in combined_df.columns if col != 'cluster']
combined_df = combined_df[cols]

# Save marker genes
combined_df.to_csv("NKT_ranked_genes_top50_by_score.csv", index=False)

#########################
# Marker genes for NKT cell annotation
genes_NKTcell = {
        'CD4+ cell': ["CD3D", "CD3E", "CD3G",'CD4'],
        'CD8+ cell': ['GZMK','CD8A','CD8B'],
        'NK cell': ["KLRD1", "NKG7", "GNLY"],
    }


# Dotplot for marker gene validation
import matplotlib.pyplot as plt
sc.pl.dotplot(adata_NKT, genes_NKTcell, groupby='leiden_NKT_res0.5', swap_axes=True, save='dotplot_NKTcell.png')
plt.show()


# Manual cell type annotation
# Map clusters to cell types
cluster_mapping = {
    '0': 'CD8+ T cells',
    '1': 'NK cells',
    '2': 'CD8+ T cells',
    '3': 'CD4+ T cells',
    '4': 'CD4+ T cells',
    '5': 'CD8+ T cells',
    '6': 'CD8+ T cells',
    '7': 'CD8+ T cells',
    '8': 'NK cells',
    '9': 'CD8+ T cells',
    '10': 'CD8+ T cells',
    '11': 'NK cells',
    '12': 'CD4+ T cells',
    '13': 'CD8+ T cells',
    '14': 'CD4+ T cells',
    '15': 'CD4+ T cells',
    '16': 'CD8+ T cells',
    '17': 'CD8+ T cells',
    '18': 'CD8+ T cells',
    '19': 'CD8+ T cells',
    '20': 'CD4+ T cells',
    '21': 'CD8+ T cells',
    '22': 'CD8+ T cells',
}


# Apply annotation
adata_NKT.obs['anno02_NKT'] = adata_NKT.obs['leiden_NKT_res0.5'].map(cluster_mapping)


# Visualize annotated cell types
sc.pl.umap(
    adata_NKT, color="anno02_NKT", legend_loc="on data", title="", frameon=False, save=".anno02_ NKT.pdf"
)

sc.pl.umap(
    adata_NKT, color="anno02_NKT", legend_loc="on data", title="", frameon=False, save=".anno02_ NKT.png"
)

sc.pl.umap(
    adata_NKT, color="anno02_NKT", legend_loc="right margin", title="", frameon=False, save=".anno02_ NKT2.png"
)


# Save annotated dataset
adata_NKT.write_h5ad('pancancer_res0.2_anno01_NKT_anno02_NKT.h5ad')


# Marker validation by dotplot
sc.pl.dotplot(adata_NKT, genes_NKTcell,'anno02_NKT', dendrogram=True,save='anno02_NKT_5.png')


########################
################### MonMac
# Extract Monocytes and Macrophages from the annotated dataset
adata_MonMac = raw_adata[raw_adata.obs['anno01']=='Monocytes and Macrophages'].copy()

# Save subset
adata_MonMac.write_h5ad('./pancancer_res0.2_anno01_Mon_Mac.h5ad')

# Load dataset
adata_MonMac=sc.read_h5ad('./pancancer_res0.2_anno01_Mon_Mac.h5ad')

# Check counts layer
adata_MonMac.layers['counts'].max()

# Normalize and log-transform
sc.pp.normalize_total(adata_MonMac, target_sum=1e4)
sc.pp.log1p(adata_MonMac)

# Store raw data
adata_MonMac.raw = adata_MonMac

# Identify highly variable genes
sc.pp.highly_variable_genes(adata_MonMac,n_top_genes=3000, flavor='seurat')

adata_MonMac.X.max()

# Scale expression values
sc.pp.scale(adata_MonMac, max_value=10)

# PCA dimensionality reduction
sc.pp.pca(adata_MonMac)

# Construct neighbor graph
sc.pp.neighbors(adata_MonMac,n_pcs = 15)

# Leiden clustering
sc.tl.leiden(adata_MonMac,flavor="igraph",n_iterations=2,resolution=0.5,key_added='leiden_MonMac_res0.5')

# UMAP visualization
#sc.tl.umap(adata_MonMac)
sc.pl.umap(adata_MonMac,color='leiden_MonMac_res0.5')

# Plot clusters on UMAP
sc.pl.umap(adata_MonMac, color=['leiden_MonMac_res0.5'], legend_loc='on data',legend_fontsize=10, title="", frameon=False)

# Save UMAP
sc.pl.umap(adata_MonMac, color=['leiden_MonMac_res0.5'], legend_fontsize=10, title="", frameon=False,save='pancancer_anno02_MonMac_celltype_0.5_1.png')

# Differential gene expression analysis
sc.tl.rank_genes_groups(adata_MonMac, 'leiden_MonMac_res0.5', method='wilcoxon')

import matplotlib.pyplot as plt

# Adjust plot margin
plt.rcParams['figure.subplot.bottom'] = 0.2

# Plot top marker genes
sc.pl.rank_genes_groups(adata_MonMac, n_genes=20, sharey=False, fontsize=10)

# Save DEG results
adata_MonMac.write("./pancancer_res0.2_anno01_counts_MonMac_res0.5_DEG.h5ad")

# Reload dataset
adata_MonMac=sc.read_h5ad("./pancancer_res0.2_anno01_counts_MonMac_res0.5_DEG.h5ad")

# Save DEG figure
ax = sc.pl.rank_genes_groups(adata_MonMac, n_genes=20, sharey=False, fontsize=10, show=False)
plt.savefig("NKT_rank_genes_groups_top20.png", dpi=300, bbox_inches='tight')


# Extract top 50 marker genes per cluster
import pandas as pd
import numpy as np

# Convert DEG results to pandas format
result = adata_MonMac.uns['rank_genes_groups']
groups = result['names'].dtype.names

ranked_genes = {}

for group in groups:
    names = result['names'][group]
    scores = result['scores'][group]

    df = pd.DataFrame({
        'gene': names,
        'score': scores
    }).sort_values('score', ascending=False).head(50)

    ranked_genes[group] = df


# Print top genes
for group, df in ranked_genes.items():
    print(f"\nGroup {group} Top 50 genes by score:")
    print(df.to_string(index=False))


# Merge all cluster markers
all_dfs = []

for group, df in ranked_genes.items():
    df_top50 = df.head(50).copy()
    df_top50["cluster"] = group
    all_dfs.append(df_top50)

combined_df = pd.concat(all_dfs, ignore_index=True)

cols = ['cluster'] + [col for col in combined_df.columns if col != 'cluster']
combined_df = combined_df[cols]

# Save marker genes
combined_df.to_csv("MonMac_ranked_genes_top50_by_score.csv", index=False)



#########################
# Marker genes for Monocyte/Macrophage annotation
genes_MonMac = {
        'Macrophages': ["CD163", "CD68",'C1QA','C1QB','C1QC','MMP9'],
        'Monocytes': ['FCN1','VCAN','APOBEC3A','THBS1']
    }

# Dotplot for marker validation
import matplotlib.pyplot as plt
sc.pl.dotplot(adata_MonMac, genes_MonMac, groupby='leiden_MonMac_res0.5', swap_axes=True, save='dotplot_MonMac_res0.5.png')
plt.show()


# Manual cell type annotation
# Map clusters to cell types
cluster_mapping = {
    '0': 'Macrophages',
    '1': 'Macrophages',
    '2': 'Macrophages',
    '3': 'Macrophages',
    '4': 'Macrophages',
    '5': 'Monocytes',
    '6': 'Macrophages',
    '7': 'Macrophages',
    '8': 'Macrophages',
    '9': 'Macrophages',
    '10': 'Macrophages',
    '11': 'Macrophages',
    '12': 'Macrophages',
    '13': 'Macrophages',
    '14': 'Macrophages',
    '15': 'Monocytes',
    '16': 'Macrophages',
    '17': 'Monocytes',
    '18': 'Monocytes',
    '19': 'Macrophages',
    '20': 'Macrophages',
    '21': 'Macrophages',
    '22': 'Macrophages',
    '23': 'Macrophages',
    '24': 'Macrophages',
}

# Apply annotation
adata_MonMac.obs['anno02_MonMac'] = adata_MonMac.obs['leiden_MonMac_res0.5'].map(cluster_mapping)

# Visualize annotated cell types
sc.pl.umap(
    adata_MonMac, color="anno02_MonMac", legend_loc="on data", title="", frameon=False, save=".anno02_MonMac.pdf"
)

sc.pl.umap(
    adata_MonMac, color="anno02_MonMac", legend_loc="on data", title="", frameon=False, save=".anno02_MonMac.png"
)

sc.pl.umap(
    adata_MonMac, color="anno02_MonMac", legend_loc="right margin", title="", frameon=False, save=".anno02_MonMac_2.png"
)

# Save annotated dataset
adata_MonMac.write_h5ad('pancancer_res0.2_anno01_MonMac_anno02_MonMac.h5ad')
