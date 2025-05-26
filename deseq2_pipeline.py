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
