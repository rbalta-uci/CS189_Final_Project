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
from preprocessing import preprocessing
if __name__ == "__main__":
    # Replace this with the path to directory where you would like results to be saved
    OUTPUT_PATH = "./output_files/synthetic_example/"
    os.makedirs(OUTPUT_PATH, exist_ok=True)  # Create path if it doesn't exist

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