import pandas as pd
import scanpy as sc
import os
from preprocessing import preprocessing

# get preprocessed data
counts_df, metadata, ann_data = preprocessing()

# create AnnData object
adata = sc.AnnData(X=counts_df.values)
adata.obs_names = counts_df.index    
adata.var_names = counts_df.columns  
adata.obs = metadata.copy()
adata.raw = adata.copy()

# set conditions from metadata
conditions = metadata["condition"].value_counts().index[:2].tolist()
condition1, condition2 = conditions[0], conditions[1]
mask = adata.obs['condition'].isin(conditions)
adata_subset = adata[mask].copy()

# differential expression analysis using Scanpy
sc.tl.rank_genes_groups(
    adata_subset,
    'condition',
    groups=[condition1],
    reference=condition2,
    method='wilcoxon',
    use_raw=True,
    key_added='wilcoxon_de'
)

# get results from Scanpy with p value and log2 fold change
wilcox_results = sc.get.rank_genes_groups_df(
    adata_subset,
    group=condition1,
    key='wilcoxon_de',
    pval_cutoff=0.05,
    log2fc_min=0.25
)

# get all results
all_wilcox_results = sc.get.rank_genes_groups_df(
    adata_subset,
    group=condition1,
    key='wilcoxon_de'
)

# format and save results
scanpy_results = all_wilcox_results.copy()
scanpy_results.columns = ['gene', 'log2FoldChange', 'pvalue', 'padj', 'stat']
scanpy_results = scanpy_results[['gene', 'log2FoldChange', 'stat', 'pvalue', 'padj']]
scanpy_results.set_index('gene', inplace=True)

output_path = "./output_files/"
os.makedirs(output_path, exist_ok=True)
scanpy_results.to_csv(os.path.join(output_path, "scanpy_results.csv"))

# compare results with DiffxPy 
try:
    diffxpy_results = pd.read_csv("output_files/examples/actual_results_diffxpy.csv")
    
    diffxpy_results_clean = diffxpy_results.copy()
    diffxpy_results_clean = diffxpy_results_clean.rename(columns={
        'pval': 'pvalue',
        'qval': 'padj',
        'log2fc': 'log2FoldChange'
    })
    diffxpy_results_clean = diffxpy_results_clean.set_index('gene')
    
    diffxpy_sig = diffxpy_results_clean[
        (diffxpy_results_clean['padj'] < 0.05) & 
        (abs(diffxpy_results_clean['log2FoldChange']) > 0.25)
    ].index
    
    scanpy_sig = wilcox_results['names'] if len(wilcox_results) > 0 else []
    overlap = set(diffxpy_sig) & set(scanpy_sig)
    
    print(f"DiffxPy significant genes: {len(diffxpy_sig)}")
    print(f"Scanpy significant genes: {len(scanpy_sig)}")
    print(f"Overlapping significant genes: {len(overlap)}")
    if len(overlap) > 0:
        print("Overlapping genes:", list(overlap)[:10])

except FileNotFoundError:
    pass
