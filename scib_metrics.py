# Import required libraries
import os
import pandas as pd
import numpy as np
import scanpy as sc

# scIB metrics benchmarking tools
from scib_metrics.benchmark import Benchmarker
from scib_metrics.benchmark import Benchmarker, BioConservation
from scib_metrics.benchmark import BatchCorrection

from matplotlib import rcParams

# Ensure text in PDF is saved as editable fonts
rcParams['pdf.fonttype'] = 42

# Load AnnData object
adata = sc.read_h5ad('data.h5ad')

# Check matrix type (dense or sparse)
type(adata.X)

# Optionally subset to HVGs only (currently commented out)
# adata = adata[:, adata.var['highly_variable']].copy()

# Set random seed for reproducibility
SEED = 42
np.random.seed(SEED)

# Subsample cells to 50,000 to reduce computational cost
sc.pp.subsample(adata, n_obs=50_000)

# Ensure HVG annotation is boolean
adata.var["highly_variable"] = np.asarray(
    adata.var["highly_variable"].astype(bool)
)

# Check available cell type labels
adata.obs['majority_voting'].unique()

# Compute PCA embedding (required for benchmarking)
sc.tl.pca(adata)

#######################################################
# Benchmark different batch correction / integration methods using scIB metrics

from scib_metrics.benchmark import Benchmarker
from scib_metrics.benchmark import Benchmarker, BioConservation
from scib_metrics.benchmark import BatchCorrection

# Optional: customize metric settings (commented out)
# biocons = BioConservation(isolated_labels=False, silhouette_label=False)
# batchcor = BatchCorrection(kbet_per_label=False)

# Initialize benchmarker
bm = Benchmarker(
    adata,
    batch_key="sample",                 # Batch annotation column
    label_key="majority_voting",        # Biological label column
    embedding_obsm_keys=[               # Embeddings to compare
        "X_pca",
        "X_harmony",
        "X_scanorama",
        "X_scVI",
        "X_BBKNN"
    ],
    n_jobs=16,                          # Parallel threads
)

# Run benchmarking
bm.benchmark()

# Retrieve results as DataFrame (no scaling)
df = bm.get_results(min_max_scale=False)

# Print results
print(df)

# Save results table
df.to_csv("batch_integration_benchmark_results.csv", index=True)

#######################################################
# Plot benchmarking results table

import matplotlib.pyplot as plt

# Generate results table plot (returns a Table object)
fig_table = bm.plot_results_table(min_max_scale=False)

# Get current figure
fig = plt.gcf()

# Save figure as high-resolution PNG and PDF
fig.savefig("benchmark_results.png",
            dpi=300, bbox_inches='tight')
fig.savefig("benchmark_results.pdf",
            dpi=300, bbox_inches='tight')

# Close figure to free memory
plt.close(fig)

#######################################################
# Save AnnData object with benchmarking results
adata.write_h5ad('./benchmark.h5ad')
