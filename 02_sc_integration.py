# Import required libraries
import omicverse as ov
# print(f"omicverse version: {ov.__version__}")
import scanpy as sc
# print(f"scanpy version: {sc.__version__}")

# Set plotting style for omicverse
ov.utils.ov_plot_set()

# Set working directory
import os
os.chdir('./data')

# Load AnnData object
adata = sc.read_h5ad('data_beforeintergration.h5ad')

# Convert expression matrix to sparse format to save memory
from scipy import sparse
adata.X = sparse.csr_matrix(adata.X)

# Save raw count matrix to layers
adata.layers["counts"] = adata.X.copy()


######## No batch correction ########

# Ensure gene names are unique
adata.var_names_make_unique()

# Basic filtering of cells and genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Backup counts layer
adata.layers["counts"] = adata.X.copy()

# Backup QC-filtered dataset
adata_bk = adata.copy()
adata.layers["counts"] = adata.X.copy()

# Normalize counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform data
sc.pp.log1p(adata)

# Store normalized data as raw
adata.raw = adata

# Identify highly variable genes (HVGs)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    batch_key="sample",
    subset=True
)

# Scale data
sc.pp.scale(adata, max_value=10)

# Perform PCA
sc.tl.pca(adata)

# Build neighbor graph
sc.pp.neighbors(adata)

# Run UMAP
sc.tl.umap(adata)

# Plot UMAP colored by batch and cell annotation
sc.pl.umap(
    adata,
    color=['batch', 'anno_Major'],
    wspace=0.5,
    save='Pcancer_anno01_unintergration.png'
)

# Save dataset
adata.write_h5ad('Pcancer_anno01_unintergration.h5ad')


##################################################
# Batch correction methods
##################################################


######## Harmony ################################

del adata2

# Copy QC-filtered dataset
adata3 = adata_bk.copy()
adata3.layers["counts"] = adata3.X.copy()

# Preprocess data using omicverse
adata3 = ov.pp.preprocess(
    adata3,
    mode='shiftlog|pearson',
    n_HVGs=3000,
    batch_key='batch'
)

# Scale data
ov.pp.scale(adata3)

# Set HVG flag
adata3.var['highly_variable'] = adata3.var['highly_variable_features'].copy()

# PCA using HVGs
sc.pp.pca(adata3, n_comps=50, mask_var="highly_variable", layer="scaled")

# Perform Harmony batch correction
sc.external.pp.harmony_integrate(adata3, 'batch')

# Build neighbor graph using Harmony embedding
sc.pp.neighbors(adata3, use_rep='X_pca_harmony', key_added='harmony_neighbours')

# Run UMAP
sc.tl.umap(adata3, neighbors_key='harmony_neighbours')

# Store UMAP embedding
adata.obsm['X_harmony'] = adata3.obsm['X_umap']

# Plot embedding
ov.utils.embedding(
    adata,
    basis='X_harmony',
    frameon='small',
    color=['batch', 'cell_type'],
    show=False
)

# Save dataset
adata.write_h5ad('Pcancer_anno01_harmony.h5ad')


######## Scanorama ################################

del adata3

adata4 = adata_bk.copy()
adata4.layers["counts"] = adata4.X.copy()

# Preprocess data
adata4 = ov.pp.preprocess(
    adata4,
    mode='shiftlog|pearson',
    n_HVGs=3000,
    batch_key='batch'
)

# Split dataset by batch
batches = adata4.obs['batch'].cat.categories.tolist()

alldata = {}
for batch in batches:
    alldata[batch] = adata4[adata4.obs['batch'] == batch,]

# Prepare datasets for Scanorama
alldata2 = dict()
for ds in alldata.keys():
    print(ds)
    alldata2[ds] = alldata[ds]

adatas = list(alldata2.values())

# Run Scanorama integration
scanorama.integrate_scanpy(adatas, dimred=50)

# Collect embeddings
scanorama_int = [ad.obsm['X_scanorama'] for ad in adatas]
all_s = np.concatenate(scanorama_int)

print(all_s.shape)

# Store Scanorama embedding
adata4.obsm["X_scanorama"] = all_s

# Build neighbors using Scanorama embedding
sc.pp.neighbors(adata4, use_rep='X_scanorama', key_added='scanorama_neighbours')

# Run UMAP
sc.tl.umap(adata4, neighbors_key='scanorama_neighbours')

# Save UMAP embedding
adata.obsm['X_scanorama'] = adata4.obsm['X_umap']

# Plot embedding
ov.utils.embedding(
    adata,
    basis='X_scanorama',
    frameon='small',
    color=['batch', 'cell_type'],
    show=False
)

# Save dataset
adata.write_h5ad('Pcancer_anno01_Scanorama.h5ad')


######## scVI ################################

del adata4

adata5 = adata_bk.copy()
adata5.layers["counts"] = adata5.X.copy()

# Preprocess data
adata5 = ov.pp.preprocess(
    adata5,
    mode='shiftlog|pearson',
    n_HVGs=3000,
    batch_key='batch'
)

# Scale data
ov.pp.scale(adata5)

# Set HVG flag
adata5.var['highly_variable'] = adata5.var['highly_variable_features'].copy()

# Setup scVI model
scvi.model.SCVI.setup_anndata(
    adata5,
    layer="counts",
    batch_key="batch"
)

# Initialize model
vae = scvi.model.SCVI(
    adata5,
    n_layers=2,
    n_latent=30,
    gene_likelihood="nb"
)

# Train model
vae.train()

# Get latent representation
adata5.obsm["X_scVI"] = vae.get_latent_representation()

# Compute MDE embedding
adata.obsm["X_mde_scVI"] = ov.utils.mde(adata5.obsm["X_scVI"])

# Plot embedding
ov.utils.embedding(
    adata,
    basis='X_mde_scVI',
    frameon='small',
    color=['batch', 'cell_type'],
    show=False
)

# Save dataset
adata.write_h5ad('Pcancer_anno01_scVI.h5ad')


######## BBKNN ################################

del adata5

adata6 = adata_bk.copy()
adata6.layers["counts"] = adata6.X.copy()

# Preprocess data
adata6 = ov.pp.preprocess(
    adata6,
    mode='shiftlog|pearson',
    n_HVGs=3000,
    batch_key='batch'
)

# Scale data
ov.pp.scale(adata6)

# Set HVG flag
adata6.var['highly_variable'] = adata6.var['highly_variable_features'].copy()

# PCA
sc.pp.pca(adata6, n_comps=50, mask_var="highly_variable", layer="scaled")

# Run BBKNN batch correction
sc.external.pp.bbknn(
    adata6,
    batch_key="batch",
    use_rep='X_pca',
    n_pcs=50,
    key_added="bbknn_neighbors"
)

# Run UMAP
sc.tl.umap(adata6, neighbors_key='bbknn_neighbors')

# Store UMAP embedding
adata.obsm['X_bbknn'] = adata6.obsm['X_umap']

# Plot embedding
ov.utils.embedding(
    adata,
    basis='X_bbknn',
    frameon='small',
    color=['batch', 'cell_type'],
    show=False
)

# Save dataset
adata.write_h5ad('Pcancer_anno01_BBKNN.h5ad')
