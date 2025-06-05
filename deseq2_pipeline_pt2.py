import os
import pickle as pkl
import scanpy as sc
import numpy as np
import pandas as pd


from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data

def get_real_data_scanpy():
    """Load real human immune cell data from scanpy"""
    # Load real human PBMC data
    adata = sc.datasets.pbmc3k_processed()
    
    # Convert to counts format for PyDESeq2
    # Get raw counts if available, otherwise use the processed data
    if 'counts' in adata.layers:
        counts_matrix = adata.layers['counts']
    else:
        # Use the main data matrix
        counts_matrix = adata.X
    
    # Convert to dense array if sparse
    if hasattr(counts_matrix, 'toarray'):
        counts_matrix = counts_matrix.toarray()
    
    # Create DataFrame (genes as rows, samples as columns)
    counts_df = pd.DataFrame(
        counts_matrix.T,  # Transpose so genes are rows
        index=adata.var.index,  # Gene names
        columns=adata.obs.index  # Cell/sample names
    )
    
    # Create metadata from the cell annotations
    metadata = adata.obs.copy()
    
    # Simplify to just a few conditions for DE analysis
    if 'louvain' in metadata.columns:
        # Use cluster assignments as conditions
        metadata['condition'] = metadata['louvain']
    elif 'leiden' in metadata.columns:
        metadata['condition'] = metadata['leiden']
    else:
        # Create artificial conditions based on first vs rest
        metadata['condition'] = ['group_A' if i < len(metadata)//2 else 'group_B' 
                                for i in range(len(metadata))]
    
    # print(f"Loaded real data: {counts_df.shape[0]} genes, {counts_df.shape[1]} cells")
    # print(f"Conditions available: {metadata['condition'].value_counts()}")
    
    return counts_df, metadata

# Use this instead of your download function
counts_df, metadata = get_real_data_scanpy()

SAVE = False  # whether to save the outputs of this notebook

if SAVE:
    # Replace this with the path to directory where you would like results to be
    # saved
    OUTPUT_PATH = "../output_files/synthetic_example"
    os.makedirs(OUTPUT_PATH, exist_ok=True)  # Create path if it doesn't exist

# counts_df = load_example_data(
#     modality="raw_counts",
#     dataset="sytnthetic",
#     debug=False,
# )

# metadata = load_example_data(
#     modality="metadata",
#     dataset="sytnthetic",
#     debug=False,
# )

print(counts_df)
print(metadata)

inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    # design="~condition",  # compare samples based on the "condition"
    # column ("B" vs "A")
    refit_cooks=True,
    inference=inference,
)

dds.fit_size_factors()

dds.obsm["size_factors"]

dds.fit_genewise_dispersions()

dds.varm["genewise_dispersions"]

dds.fit_dispersion_trend()
dds.uns["trend_coeffs"]
dds.varm["fitted_dispersions"]

dds.fit_dispersion_prior()
# print(
#     f"logres_prior={dds.uns['_squared_logres']}, sigma_prior={dds.uns['prior_disp_var']}"
# )

dds.fit_MAP_dispersions()
dds.varm["MAP_dispersions"]
dds.varm["dispersions"]

dds.fit_LFC()
dds.varm["LFC"]

dds.refit()

ds = DeseqStats(
    dds,
    contrast=np.array([0, 1]),
    #contrast=["condition", "B", "A"],
    alpha=0.05,
    cooks_filter=True,
    independent_filter=True,
)

ds.run_wald_test()
ds.p_values

if ds.independent_filter:
    ds._independent_filtering()
else:
    ds._p_value_adjustment()

ds.padj

ds.summary()

if SAVE:
    with open(os.path.join(OUTPUT_PATH, "stat_results_detailed_pipe.pkl"), "wb") as f:
        pkl.dump(ds, f)

ds.lfc_shrink(coeff="condition_B_vs_A")

if SAVE:
    with open(os.path.join(OUTPUT_PATH, "shrunk_results_detailed_pipe.pkl"), "wb") as f:
        pkl.dump(ds, f)