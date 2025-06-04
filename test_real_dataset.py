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
    

    ann_data.obs["cell_type"] = ann_data.obs["Factor Value[inferred cell type - authors labels]"].astype("category")
    
    mask = ~ann_data.obs["cell_type"].isna() #get rid of nan columns
    ann_data = ann_data[mask].copy()

    
    ann_data.X = ann_data.X.toarray() #make dense

    ann_data = ann_data[:1000, :]     # First 1000 cells

    min_cells_per_type = 10
    valid_cell_types = ann_data.obs["cell_type"].value_counts()
    valid_cell_types = valid_cell_types[valid_cell_types >= min_cells_per_type].index

    mask = ann_data.obs["cell_type"].isin(valid_cell_types)
    ann_data = ann_data[mask].copy()


    ann_data.obs["cell_type"] = ann_data.obs["cell_type"].cat.remove_unused_categories()

    sc.pp.filter_genes(ann_data, min_cells=5)

    print("Filtered shape:", ann_data.shape)
    print("Genes with all zeros:", (ann_data.X.sum(axis=0) == 0).sum())

    full_formula_loc = "~ 1 + cell_type"

    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    ann_data.X = np.nan_to_num(ann_data.X, nan=0.0, posinf=0.0, neginf=0.0)

    nonzero_genes = (ann_data.X.sum(axis=0) > 0).A1 if hasattr(ann_data.X, "A1") else (ann_data.X.sum(axis=0) > 0)
    ann_data = ann_data[:, nonzero_genes]
    ann_data.X = np.asarray(ann_data.X, dtype=np.float32)
    
    # Convert to DataFrame for easier group-wise variance calculation
    X_df = pd.DataFrame(ann_data.X, columns=ann_data.var_names)
    X_df["cell_type"] = ann_data.obs["cell_type"].values

    # Compute variance per gene within each group
    group_variance = X_df.groupby("cell_type").var()

    # Identify genes with nonzero variance across all groups
    nonzero_var_genes = group_variance.columns[(group_variance > 1e-6).all(axis=0)]

    # Filter ann_data to include only those genes
    ann_data = ann_data[:, nonzero_var_genes]

    test = de.test.wald(
        data=ann_data,
        formula_loc="~ 1 + cell_type",
        factor_loc_totest="cell_type",
        sample_description=ann_data.obs,
        gene_names=ann_data.var_names.tolist()
    )

    print(test.summary())
    test = test.summary()
    print("Number of genes tested:", test.shape[0])
    print("Significant genes:", (test['pval'] < 0.05).sum())
    plt.scatter(x=test['log2fc'],y=test['pval'].apply(lambda x:-np.log10(x)),s=1)
    plt.show()