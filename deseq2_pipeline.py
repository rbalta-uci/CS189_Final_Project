# Steps to conduct differential analysis

# implement our libraries 

# upload our data into our envirionnment
# this includes our raw data and our metadata

import os
import pickle as pkl

import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# Replace this with the path to directory where you would like results to be saved
OUTPUT_PATH = "../output_files/synthetic_example/"
os.makedirs(OUTPUT_PATH, exist_ok=True)  # Create path if it doesn't exist

# Replace this with the path to your dataset
DATA_PATH = "https://raw.githubusercontent.com/owkin/PyDESeq2/main/datasets/synthetic/"

counts_df = pd.read_csv(os.path.join(DATA_PATH, "test_counts.csv"), index_col=0)
print(counts_df)

counts_df = counts_df.T
print(counts_df)

metadata = pd.read_csv(os.path.join(DATA_PATH, "test_metadata.csv"), index_col=0)
print(metadata)

samples_to_keep = ~metadata.condition.isna()
counts_df = counts_df.loc[samples_to_keep]
metadata = metadata.loc[samples_to_keep]

genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
counts_df = counts_df[genes_to_keep]

inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    refit_cooks=True,
    inference=inference,
)

dds.deseq2()

print(dds)

with open(os.path.join(OUTPUT_PATH, "result_dds.pkl"), "wb") as f:
    pkl.dump(dds, f)
    
print(dds.varm["dispersions"])

print(dds.varm["LFC"])

ds = DeseqStats(dds, contrast=["condition", "B", "A"], inference=inference)

ds.summary()

ds.results_df.to_csv(os.path.join(OUTPUT_PATH, "results.csv"))