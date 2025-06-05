import pandas as pd
import numpy as np
import scanpy as sc
def preprocessing():
    ann_data = sc.datasets.ebi_expression_atlas("E-MTAB-9543")

    ann_data.obs["cell_type"] = ann_data.obs["Factor Value[inferred cell type - authors labels]"].astype("category")

    # Drop NA cell types
    ann_data = ann_data[~ann_data.obs["cell_type"].isna()].copy()
    sc.pp.filter_genes(ann_data, min_cells=10)
    sc.pp.filter_cells(ann_data, min_counts=1000)
    # Convert to DataFrame
    counts_df = pd.DataFrame.sparse.from_spmatrix( ann_data.X.astype(np.float32), index=ann_data.obs_names, columns=ann_data.var_names)
    counts_df = counts_df.astype(np.int32)
    # Metadata
    metadata = ann_data.obs.copy()
    metadata["condition"] = metadata["cell_type"]  
    return counts_df, metadata, ann_data

def relevant_genes(file, output):
    with open(file, 'r') as f:
        with open(output, 'w') as r:
            r.write(f.readline())
            for line in f:
                line_split = line.split(',')
                
                try:
                    if float(line_split[5])< .05:
                        r.write(line)
                except:
                    continue

if __name__ == "__main__":
    relevant_genes("./output_files/synthetic_example/results.csv", "./output_files/actual_results.csv")
