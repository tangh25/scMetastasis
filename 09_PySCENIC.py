# PySCENIC workflow for single-cell regulatory network analysis

import os
import math
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss

# ----------------- Settings -----------------
warnings.simplefilter(action="ignore", category=Warning)
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80)

projpath = './analysis'
chapter = "pySCENIC"

# Create output folder
os.chdir(projpath)
if not os.path.exists(chapter):
    os.makedirs(chapter)
os.chdir(chapter)

# ----------------- Load data -----------------
# Load single-cell AnnData object
adata = sc.read_h5ad(
    "./adata.h5ad"
)

# Optional: remove unwanted cell types (commented out)
# remove_types = ["Mesenchymal stem-like cells", "Pericyte", "SMC"]
# adata = adata[~adata.obs["anno_Sub"].isin(remove_types)].copy()

print(f"Data shape: {adata.shape}")
print(f"Max expression: {adata.X.max()}")

# Convert AnnData to DataFrame: cells x genes
ex_matrix = adata.to_df()

# ----------------- Load TFs & motif databases -----------------
tf_names = load_tf_names(os.path.join(projpath, "pySCENIC/pyscenic/allTFs_hg38.txt"))

hg38_10kbp = os.path.join(projpath,
                          "pySCENIC/pyscenic/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
hg38_500bp = os.path.join(projpath,
                          "pySCENIC/pyscenic/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")

dbs = [RankingDatabase(fname=hg38_10kbp, name="hg38_10kbp"),
       RankingDatabase(fname=hg38_500bp, name="hg38_500bp")]

# TF motif annotation
tf_anno = os.path.join(projpath, "pySCENIC/pyscenic/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")

# ----------------- Co-expression network -----------------
# Compute adjacency matrix (can take hours for large datasets)
adjacencies = grnboost2(
    ex_matrix.iloc[:200, :],  # subset for testing; replace with ex_matrix for full dataset
    tf_names=tf_names,
    client_or_address='local',
    verbose=True
)
adjacencies.to_csv("adjacencies.csv", index=False)

# ----------------- Build modules -----------------
adjacencies = pd.read_csv("adjacencies.csv")
modules = list(modules_from_adjacencies(adjacencies, ex_matrix))

# ----------------- Motif pruning -----------------
with ProgressBar():
    df = prune2df(dbs, modules, tf_anno)
df.to_csv("motifs.csv")
df.head()

# ----------------- Compute regulons -----------------
regulons = df2regulons(df)
# Remove regulons with few target genes
regulons = [r for r in regulons if len(r.gene2weight) >= 10]

# Compute AUCell scores (regulon activity per cell)
auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
auc_mtx.to_csv("auc_mtx.csv")

# ----------------- Binarize regulon activity -----------------
auc_mtx = pd.read_csv("auc_mtx.csv", index_col=0)
bin_mtx, thresholds = binarize(auc_mtx, num_workers=8)
bin_mtx.to_csv("bin_mtx.csv")
thresholds.to_frame().rename(columns={0: 'threshold'}).to_csv("bin_thresholds.csv")

# ----------------- Regulon specificity scores -----------------
rss_celltype = regulon_specificity_scores(auc_mtx, adata.obs['anno_Sub'])
rss_celltype.to_csv("rss_celltype.csv")
rss_celltype

# ----------------- Plot RSS per cell type -----------------
cats = sorted(adata.obs['anno_Sub'].unique())
ncols = 5
nrows = math.ceil(len(cats) / ncols)

fig = plt.figure(figsize=(4*ncols, 3*nrows))
for idx, c in enumerate(cats):
    ax = fig.add_subplot(nrows, ncols, idx + 1)
    plot_rss(rss_celltype, c, top_n=5, max_n=None, ax=ax)
    x = rss_celltype.T[c]
    dy = (x.max() - x.min()) * 0.05
    ax.set_ylim(x.min() - dy, x.max() + dy)
    ax.set_title(c, fontsize=14)
    ax.set_xlabel('')
    ax.set_ylabel('')
    for t in ax.texts:
        t.set_fontsize(7)
    adjust_text(ax.texts, only_move='x+y', arrowprops=dict(arrowstyle='-', color='lightgrey'))

# Shared axis labels
fig.text(0.5, 0.02, 'Regulon', ha='center', va='center', size='x-large')
fig.text(0.06, 0.5, 'Regulon specificity score (RSS)', ha='center',
         va='center', rotation='vertical', size='x-large')

plt.tight_layout(rect=[0.1, 0.05, 1, 0.95])
plt.savefig("rss_plot3.pdf", format="pdf", dpi=300)
plt.show()

# ----------------- Heatmaps -----------------
df = pd.read_csv("./MAC/rss_celltype.csv", index_col=0)
df.columns = [x.replace("(+)", "") for x in df.columns]


# Full heatmap
plt.figure(figsize=(20, 6))
sns.heatmap(df, cmap="YlGnBu", linewidths=0.3, cbar_kws={"label": "RSS value"})
plt.title("Regulon Specificity Scores (RSS) across Macrophage Subtypes", fontsize=16)
plt.xlabel("Transcription Factors", fontsize=12)
plt.ylabel("Cell Types", fontsize=12)
plt.tight_layout()
plt.savefig("rss_heatmap_all.png", dpi=300)
plt.close()

# Top 50 TFs by mean RSS
top_tfs = df.mean(axis=0).sort_values(ascending=False).head(50).index
df_top = df[top_tfs]
plt.figure(figsize=(16, 6))
sns.heatmap(df_top, cmap="coolwarm", linewidths=0.3, cbar_kws={"label": "RSS value"})
plt.title("Top 50 TFs by mean RSS", fontsize=16)
plt.xlabel("Transcription Factors", fontsize=12)
plt.ylabel("Cell Types", fontsize=12)
plt.tight_layout()
plt.savefig("rss_heatmap_top50.png", dpi=300)
plt.close()

# Clustered heatmap (rows and columns)
sns.clustermap(df, cmap="vlag", standard_scale=1,
               figsize=(14, 8),
               linewidths=0.2,
               cbar_kws={"label": "Scaled RSS"})
plt.savefig("rss_clustermap.png", dpi=300)
plt.close()

print("All plots completed!")
