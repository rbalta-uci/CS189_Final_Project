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

def run_diffxpy():
    count_df, metadata, anndata = preprocessing()
    if not isinstance(anndata.X, np.ndarray):
        anndata.X = anndata.X.toarray()

    test = de.test.wald(
        data=anndata,
        formula_loc="~ 1 + cell_type",
        factor_loc_totest="cell_type",  
        sample_description=metadata,
        gene_names = anndata.var_names
    )
    test = test.summary()
    test.to_csv("./output_files/diffxpy.csv",index=False)

def create_graph_and_filter_data(results):
    df = pd.read_csv(results)
    df = df[df["pval"].notna and df["pval"]<.05]
    df.to_csv("./output_files/actual_results_diffxpy.csv",index=False)
    plt.scatter(x=np.log2(df['mean']), y=df['log2fc'])
    plt.xlabel("Mean Expression", fontsize=12)
    plt.ylabel("Log2 Fold Change", fontsize=12)
    plt.title("MA Plot of Significant Genes", fontsize=14)
    plt.savefig('./output_files/plot_diffxpy.png')
    plt.show()


if __name__ == "__main__":    
    #run_diffxpy()
    create_graph_and_filter_data("./output_files/diffxpy.csv")
