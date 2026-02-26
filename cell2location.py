
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
import sys
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
from matplotlib import rcParams

rcParams['pdf.fonttype'] = 42  # enables correct plotting of text for PDFs
import os

data_in_path = "cell2location"
os.chdir(data_in_path)

results_folder = './results/analysis_major/'
# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

# adata_vis = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
adata_vis = sc.read_visium(path="./", count_file='./filtered_feature_bc_matrix.h5', library_id="GSM8651955"
                           , load_images=True, source_image_path="./spatial/")
adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]
adata_vis.var['SYMBOL'] = adata_vis.var_names
# adata_vis.var.set_index('gene_ids', drop=True, inplace=True)
sc.pl.spatial(adata=adata_vis, color='PTPRC', gene_symbols='SYMBOL')
sc.pl.spatial(adata=adata_vis, color='PTPRC', gene_symbols='SYMBOL', save='PTPRC.png')
import squidpy as sq

adata_vis.var['SYMBOL'] = adata_vis.var.index
# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]
# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

adata_vis.write_h5ad('adata_vis.h5ad')

adata_ref = sc.read_h5ad("data.h5ad")
del adata_ref.raw

adata_ref.var['SYMBOL'] = adata_ref.var.index
adata_ref.X = adata_ref.layers['counts'].copy()
adata_ref.X.max()
adata_ref.var_names_make_unique()

sc.settings.set_figure_params(dpi=100, dpi_save=300, figsize=(4, 4))
sc.pl.umap(adata_ref, color='Sub_celltype')

adata_ref.var['Mt'] = [gene.startswith('MT-') for gene in adata_ref.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_ref.obsm['Mt'] = adata_ref[:, adata_ref.var['Mt'].values].X.toarray()  
adata_ref = adata_ref[:, ~adata_ref.var['Mt'].values]

from cell2location.utils.filtering import filter_genes

selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
# filter the object
adata_ref = adata_ref[:, selected].copy()


import numpy as np
np.random.seed(42)
cell_types = adata_ref.obs['Sub_celltype'].unique()
selected_indices = []
for ct in cell_types:
    indices = adata_ref.obs[adata_ref.obs['Sub_celltype'] == ct].index
    if len(indices) <= 300:
        selected = indices
    else:
        selected = np.random.choice(indices, size=300, replace=False)
    selected_indices.extend(selected)
adata_ref = adata_ref[selected_indices].copy()


# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                                                   # 10X reaction / sample / batch
                                                   batch_key='sample',
                                                   # cell type, covariate used for constructing signatures
                                                   labels_key='Sub_celltype',  # cell_type
                                                   # multiplicative technical effects (platform, 3' vs 5', donor effect)
                                                   categorical_covariate_keys=None 
                                                   )

# create the regression model
from cell2location.models import RegressionModel

mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()

adata_vis.X = np.round(adata_vis.X).astype(int)
adata_ref.X = np.round(adata_ref.X).astype(int)

mod.train(max_epochs=1000, batch_size=2500, accelerator='gpu')

history_df = mod.history_['elbo_train']

print(history_df.head()) 
print(history_df.columns)  
import matplotlib.pyplot as plt
plt.figure(figsize=(8, 5))
plt.plot(history_df.index, history_df["elbo_train"], label="ELBO (train)")
plt.xlabel("Epoch")
plt.ylabel("ELBO")
plt.title("Training history of RegressionModel")
plt.legend()
plt.grid(True)
plt.savefig("01-mod_training_elbo_GSM8651955.png", dpi=300, bbox_inches="tight")
plt.show()

#import matplotlib.pyplot as plt
mod.plot_history(20)
plt.savefig("01-mod.plot_history.png", format='png', bbox_inches='tight')  
plt.show()  

#In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}  # , 'use_cpu': True
)

# mod.export_posterior(adata_ref)
# Save model
mod.save(f"{ref_run_name}", overwrite=True)
# Save anndata object with results
adata_file = f"{ref_run_name}/sc_GSM8651955.h5ad"
adata_ref.write(adata_file)
adata_file

# Save model
mod.save("./", overwrite=True)
mod.plot_QC()
plt.savefig('02-mod.plot_QC.png') 

# Save anndata object with results
adata_ref.write("./sc_cell2location_GSM8651955.h5ad")
# compute the 5%, 50% and 95% quantiles of the posterior distribution directly rather than using 1000 samples from the distribution (or any other quantiles)
adata_ref = mod.export_posterior(
    adata_ref, use_quantiles=True,
    # choose quantiles
    add_to_varm=["q05", "q50", "q95", "q0001"],
    sample_kwargs={'batch_size': 2500}  # , 'use_gpu': True
)


mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                                          for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                              for i in adata_ref.uns['mod']['factor_names']]].copy()

inf_aver.columns = adata_ref.uns['mod']['factor_names']
print(inf_aver.iloc[0:12, 0:12])

####Cell2location: spatial mapping
# find shared genes and subset both anndata and reference signatures
adata_vis.var_names_make_unique()
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis,
                                                 batch_key="sample") 
print(adata_vis.X[:10, :10])
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    N_cells_per_location=30,
    detection_alpha=20
)
mod.view_anndata_setup()

mod.train(max_epochs=10000,  # 30000
          batch_size=None,
          train_size=1,
          accelerator='gpu'
          )
mod.plot_history()

history_df = mod.history_['elbo_train']

print(history_df.head())   
print(history_df.columns) 
import matplotlib.pyplot as plt
plt.figure(figsize=(8, 5))
plt.plot(history_df.index, history_df["elbo_train"], label="ELBO (train)")
plt.xlabel("Epoch")
plt.ylabel("ELBO")
plt.title("Training history of RegressionModel")
plt.legend()
plt.grid(True)
plt.savefig("03-mod_training_elbo_GSM8651955.png", dpi=300, bbox_inches="tight")
plt.show()

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);
plt.savefig('mod.plot_history2.png')  
plt.show()  
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}  # , 'accelerator': 'gpu'
)
# Save model
run_name = './'
mod.save(f"{run_name}", overwrite=True)

adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

mod.plot_QC()
plt.savefig('04-mod.plot_QC_GSM8651955.png') 
plt.show()

adata_vis.obsm
adata_vis.obsm['means_cell_abundance_w_sf']
import pandas as pd

pd.DataFrame(adata_vis.obsm['q05_cell_abundance_w_sf']).to_csv("./st_cell2location_GSM8651955_res.csv")
adata_vis.uns['spatial'].keys()  

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf'] 

adata_vis.write_h5ad('./adata_vis_GSM8651955.h5ad')
adata_ref.write_h5ad('./adata_ref_GSM8651955.h5ad')
adata_vis=sc.read_h5ad('./adata_vis_GSM8651955.h5ad')

from cell2location.utils import select_slide

slide = select_slide(adata_vis, 'GSM8651955')

# Now we use cell2location plotter that allows showing multiple cell types in one panel
from cell2location.plt import plot_spatial


clust_labels = ['Macro_LYVE1','Fib_HOPX']   
clust_col = ['' + str(i) for i in clust_labels]  # in case column names differ from labels

slide = select_slide(adata_vis, 'GSM8651955')

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
        show_img=True,
        # 'fast' (white background) or 'dark_background'
        style='fast',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=6,
        colorbar_position='right'
    )
plt.savefig('sc.pl.spatial_Macro_LYVE1_HOPX+ Fib.pdf')  