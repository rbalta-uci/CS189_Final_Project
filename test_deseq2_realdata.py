import os
import pickle as pkl
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import diffxpy.api as de
import matplotlib.pylab as plt
import numpy as np
import scanpy as sc
def preprocessing():
    os.makedirs(OUTPUT_PATH, exist_ok=True)  # Create path if it doesn't exist

    ann_data = sc.datasets.ebi_expression_atlas("E-MTAB-9543")

    ann_data.obs["cell_type"] = ann_data.obs["Factor Value[inferred cell type - authors labels]"].astype("category")

    # Drop NA cell types
    ann_data = ann_data[~ann_data.obs["cell_type"].isna()].copy()
    sc.pp.filter_genes(ann_data, min_cells=10)
    sc.pp.filter_cells(ann_data, min_counts=1000)  # optional
    # Convert to DataFrame
    counts_df = pd.DataFrame.sparse.from_spmatrix( ann_data.X.astype(np.float32), index=ann_data.obs_names, columns=ann_data.var_names)
    counts_df = counts_df.astype(np.int32)
    # Metadata
    metadata = ann_data.obs.copy()
    metadata["condition"] = metadata["cell_type"]  # Add 'condition' column
    return counts_df, metadata
if __name__ == "__main__":
    # Replace this with the path to directory where you would like results to be saved
    OUTPUT_PATH = "./output_files/synthetic_example/"

    counts_df, metadata = preprocessing()
    # DESeq2 Setup
    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        refit_cooks=True,
        inference=inference,
    )

    dds.deseq2()

    # Contrast setup (automatically pick top 2 conditions)
    conditions = metadata["condition"].value_counts().index[:2].tolist()
    print("Using contrast between:", conditions)

    ds = DeseqStats(dds, contrast=["condition", conditions[0], conditions[1]], inference=inference)
    ds.summary()

    ds.results_df.to_csv(os.path.join(OUTPUT_PATH, "results.csv"))