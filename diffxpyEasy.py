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
import diffxpy.api as de
if __name__ == "__main__":
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


    full_formula_loc = '~ 1 + condition'
    reduced_formula_loc = '~ 1'
    sample_description=metadata

    print(type(metadata))
    print("Type of full_formula_loc:", type(full_formula_loc))
    print("Type of reduced_formula_loc:", type(reduced_formula_loc))
    print(metadata.columns.tolist())
    test = de.test.wald(counts_df.values, 
                    formula_loc=full_formula_loc, 
                    sample_description=metadata,
                    factor_loc_totest = "condition",
                    gene_names= counts_df.columns.tolist())

    print(test.summary().head())
