### Basic workflows
### Preprocessing and clustering

# Import required libraries
import scanpy as sc
import numpy as np
import anndata as ad
import pandas as pd
import os
import matplotlib.pyplot as plt

# Set verbosity level for scanpy
sc.settings.verbosity = 3
sc.logging.print_header()

# Set global plotting parameters (white background, 300 dpi)
sc.settings.set_figure_params(dpi=300, facecolor="white")

# Set input data path
data_in_path = "./"

# Change working directory
os.chdir(data_in_path)
print(os.getcwd())

# Read h5ad file
file_path = './'
adata = sc.read_h5ad(file_path)
print("Data loaded successfully")

# Get file name without extension
file_name = os.path.basename(file_path).split('_')[0]

# Create output folder
folder_path = os.path.join(data_in_path, file_name)
os.makedirs(folder_path, exist_ok=True)
print(f"Create folder: {folder_path}")

# Check data matrix shape
adata.to_df().shape

# Display first 5 rows and 5 columns of expression matrix
adata.to_df().iloc[0:5,0:5]

# Preview metadata
adata.obs.head()

# Print unique sample IDs
print(adata.obs['orig.ident'].unique())

# Check unique values in 'orig.ident'
print(adata.obs['orig.ident'].unique())

# Confirm unique values again
print(adata.obs['orig.ident'].unique())

# Print all unique raw cell types
print(list(adata.obs['raw_celltype'].unique()))

# Add sample grouping information
adata.obs['sample'] = adata.obs['orig.ident']

# Extract sample name by splitting 'orig.ident' with "_"
# Map the first part to sample group using predefined mapping
adata.obs['sample'] = adata.obs['orig.ident'].str.split('_').str[0].map(mapping)

# Display first few rows to verify new 'sample' column
print(adata.obs[['orig.ident', 'sample']].head())

# Check unique sample groups
adata.obs['sample'].unique()


# =========================
# Quality Control (QC)
# =========================

# Identify mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Identify hemoglobin genes
adata.var['hb'] = adata.var_names.str.contains(("^HB[^(P)]"))

# Calculate QC metrics (mitochondrial and hemoglobin gene percentage)
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','hb'], percent_top=None, log1p=False, inplace=True)

# Visualize QC metrics for each sample
sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_hb'],
    jitter=0.4,
    groupby='sample',
    rotation=45,
    save="_before_qc_p1.pdf"
)

# =========================
# Filtering
# =========================

# Plot top 20 highly expressed genes
sc.pl.highest_expr_genes(adata, n_top=20, show=False)
plt.savefig('./pre_qc_total_counts_top20_genes.pdf', bbox_inches='tight')

# Filter cells and genes based on QC thresholds
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Remove cells with high mitochondrial or hemoglobin gene percentage
adata = adata[adata.obs['pct_counts_mt'] <= 20, :]
adata = adata[adata.obs['pct_counts_hb'] <= 1, :]

print(adata.n_obs, adata.n_vars)

# =========================
# Doublet detection using Scrublet
# =========================

# Save raw data before doublet filtering
adata.raw = adata

# Import Scrublet
import scrublet as scr

# Run doublet detection
scrub = scr.Scrublet(adata.raw.X, expected_doublet_rate=0.06)
result = scrub.scrub_doublets(verbose=False, n_prin_comps=20)

# Plot doublet score histogram
scrub.plot_histogram()

# Add doublet prediction results to metadata
adata.obs['doublet_scores'] = result[0]
adata.obs['predicted_doublets'] = result[1]

# Count predicted doublets
adata.obs.predicted_doublets.value_counts()

# Convert boolean to string for plotting
adata.obs['predicted_doublets'] = adata.obs["predicted_doublets"].astype(str)

# Visualize gene number distribution for doublets vs singlets
sc.pl.violin(
    adata,
    'n_genes_by_counts',
    jitter=0.4,
    groupby='predicted_doublets',
    show=False
)
plt.savefig('./predicted_doublets_violin_plot.pdf', bbox_inches='tight')


# =========================
# Feature selection and dimensionality reduction
# =========================

# Identify highly variable genes (HVGs)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3", subset=False)

# Subset to HVGs
adata = adata[:, adata.var.highly_variable]

# Regress out unwanted effects (library size and mitochondrial content)
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# Scale data
sc.pp.scale(adata, max_value=10)

# Perform PCA
sc.tl.pca(adata, svd_solver='arpack')

# Compute neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Run UMAP
sc.tl.umap(adata)


# =========================
# Visualization
# =========================

# Plot UMAP colored by doublet scores and predictions
lis = ['doublet_scores','predicted_doublets']
for i in lis:
    sc.pl.umap(adata, color=i, show=False)
    plt.savefig('./' + i + '_doublet_scores.pdf', bbox_inches='tight')


# =========================
# Remove doublets
# =========================

# Restore raw data
adata = adata.raw.to_adata()

# Keep only singlets
adata = adata[adata.obs['predicted_doublets'] == 'False', :]

print(adata.shape)


# =========================
# QC after filtering
# =========================

sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_hb'],
    jitter=0.4,
    groupby='orig.ident',
    rotation=45,
    save="_after_qc_p2.pdf"
)

# Save processed data
adata.write_h5ad('./afterqc.h5ad', compression='gzip')
