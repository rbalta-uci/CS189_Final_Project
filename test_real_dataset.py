import os
import pickle as pkl
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import diffxpy.api as de
import matplotlib.pylab as plt
import numpy as np
import scipy.io
import scanpy as sc
if __name__ == "__main__":
    # Replace this with the path to directory where you would like results to be saved
    #OUTPUT_PATH = "../output_files/synthetic_example/"
    #os.makedirs(OUTPUT_PATH, exist_ok=True)  # Create path if it doesn't exist

    # Replace this with the path to your dataset
    #DATA_PATH = "https://raw.githubusercontent.com/owkin/PyDESeq2/main/datasets/synthetic/"

    #counts_df = pd.read_csv(os.path.join(DATA_PATH, "test_counts.csv"), index_col=0)

    ann_data = sc.datasets.ebi_expression_atlas("E-MTAB-9543")
    counts_df = ann_data.to_df()
    metadata =pd.read_csv('ExpDesign-E-MTAB-9543.tsv', sep='\t')
      
    print(ann_data.shape)
    #print(metadata.columns)
   # Filter metadata samples with non-NA cell type label
    samples_to_keep = counts_df['Factor Value[inferred cell type - authors labels]'].notna()
    filtered_sample_ids = counts_df.index[samples_to_keep]

    # Subset counts_df and metadata by matching sample IDs explicitly:
    counts_df = counts_df.loc[filtered_sample_ids, :]  # rows=samples
    metadata = counts_df.loc[filtered_sample_ids, :]

    # Filter genes with at least one count across filtered samples
    genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 1]
    counts_df = counts_df[genes_to_keep]

    # Check shapes again
    print("Counts shape:", counts_df.shape)    # (samples, genes)
    print("Metadata shape:", metadata.shape)   # (samples, metadata columns)

    full_formula_loc = "~ 1 + 'Factor Value[inferred cell type - authors labels]'"
    test = de.test.wald(counts_df.values, 
                    formula_loc=full_formula_loc, 
                    sample_description=metadata,
                    factor_loc_totest = "Factor Value[inferred cell type - authors labels]",
                    gene_names= counts_df.columns.tolist())

    print(test.summary())
    test = test.summary()
    print("Number of genes tested:", test.shape[0])
    print("Significant genes:", (test['pval'] < 0.05).sum())
    plt.scatter(x=test['log2fc'],y=test['pval'].apply(lambda x:-np.log10(x)),s=1)
    plt.show()