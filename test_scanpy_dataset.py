import pandas as pd
import scanpy as sc
import os
from preprocessing import preprocessing  

counts_df, metadata = preprocessing()

# get cell and gene names from data
index_sample = str(counts_df.index[0]) if len(counts_df.index) > 0 else ""
columns_sample = str(counts_df.columns[0]) if len(counts_df.columns) > 0 else ""

if "ENSG" in index_sample or ("-" in index_sample and len(index_sample) > 15):
    count_matrix = counts_df
    cell_names = counts_df.index
    gene_names = counts_df.columns
elif "ENSG" in columns_sample or ("-" in columns_sample and len(columns_sample) > 15):
    count_matrix = counts_df.T
    cell_names = counts_df.columns
    gene_names = counts_df.index
else:
    if counts_df.shape[0] < counts_df.shape[1]:
        count_matrix = counts_df.T
        cell_names = counts_df.columns
        gene_names = counts_df.index
    else:
        count_matrix = counts_df
        cell_names = counts_df.index
        gene_names = counts_df.columns

# create AnnData object
adata = sc.AnnData(X=count_matrix.values)
adata.obs_names = cell_names
adata.var_names = gene_names
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
import os
os.makedirs(output_path, exist_ok=True)
scanpy_results.to_csv(os.path.join(output_path, "scanpy_results.csv"))

# compare results with DESeq2
try:
    deseq2_results = pd.read_csv("./output_files/results.csv", index_col=0)
    deseq2_sig = deseq2_results[deseq2_results['padj'] < 0.05].index
    scanpy_sig = wilcox_results['names'] if len(wilcox_results) > 0 else []
    overlap = set(deseq2_sig) & set(scanpy_sig)
    
    print(f"DESeq2 significant genes: {len(deseq2_sig)}")
    print(f"Scanpy significant genes: {len(scanpy_sig)}")
    print(f"Overlapping significant genes: {len(overlap)}")
    if len(overlap) > 0:
        print("Overlapping genes:", list(overlap)[:10])

except FileNotFoundError:
    pass