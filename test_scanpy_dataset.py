import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from matplotlib_venn import venn2
from preprocessing import preprocessing  # Import your preprocessing function

warnings.simplefilter(action="ignore", category=Warning)
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=80)

print("="*50)
print("LOADING DATA WITH PREPROCESSING")
print("="*50)

# Use your preprocessing function (same as DESeq2 analysis)
counts_df, metadata = preprocessing()

print(f"Counts shape: {counts_df.shape}")
print(f"Metadata shape: {metadata.shape}")
print(f"Metadata columns: {metadata.columns.tolist()}")

# Check the data
print("\nFirst few rows of counts:")
print(counts_df.head())
print("\nFirst few rows of metadata:")
print(metadata.head())
print("\nCondition distribution:")
print(metadata["condition"].value_counts())

# Convert to AnnData object for scanpy
print("\n" + "="*50)
print("CONVERTING TO ANNDATA")
print("="*50)

# Check dimensions and alignment
print(f"Counts DataFrame shape: {counts_df.shape}")
print(f"Metadata shape: {metadata.shape}")

# Determine correct orientation by looking at the actual names
# Check if index looks like cell names (contain dashes/barcodes) or gene names (ENSG)
index_sample = str(counts_df.index[0]) if len(counts_df.index) > 0 else ""
columns_sample = str(counts_df.columns[0]) if len(counts_df.columns) > 0 else ""

print(f"Index sample: {index_sample}")
print(f"Columns sample: {columns_sample}")

# Cell barcodes typically contain dashes, gene IDs typically start with ENSG
if "ENSG" in index_sample or ("-" in index_sample and len(index_sample) > 15):
    # Index contains cell names, columns contain gene names
    print("Cells are rows, genes are columns - correct orientation")
    count_matrix = counts_df
    cell_names = counts_df.index
    gene_names = counts_df.columns
elif "ENSG" in columns_sample or ("-" in columns_sample and len(columns_sample) > 15):
    # Columns contain cell names, index contains gene names
    print("Genes are rows, cells are columns - transposing")
    count_matrix = counts_df.T
    cell_names = counts_df.columns
    gene_names = counts_df.index
else:
    # Fall back to size comparison
    if counts_df.shape[0] < counts_df.shape[1]:
        print("Based on size: genes are rows, cells are columns - transposing")
        count_matrix = counts_df.T
        cell_names = counts_df.columns
        gene_names = counts_df.index
    else:
        print("Based on size: cells are rows, genes are columns - no transpose needed")
        count_matrix = counts_df
        cell_names = counts_df.index
        gene_names = counts_df.columns

print(f"Final count matrix shape: {count_matrix.shape}")
print(f"Number of cells: {len(cell_names)}")
print(f"Number of genes: {len(gene_names)}")

# Ensure metadata matches cell names
if len(metadata) != len(cell_names):
    print(f"WARNING: Metadata length ({len(metadata)}) doesn't match cell count ({len(cell_names)})")
    # Try to align by index
    if set(metadata.index) == set(cell_names):
        print("Reordering metadata to match cell order")
        metadata = metadata.loc[cell_names]
    else:
        print("ERROR: Cannot align metadata with cells")
        print(f"Metadata index sample: {metadata.index[:5].tolist()}")
        print(f"Cell names sample: {cell_names[:5].tolist()}")
        raise ValueError("Metadata and count matrix don't align")

# Create AnnData object
adata = sc.AnnData(X=count_matrix.values)
adata.obs_names = cell_names
adata.var_names = gene_names

# Add metadata to observations
adata.obs = metadata.copy()

print(f"AnnData shape: {adata.shape}")
print(f"Observations (cells/samples): {adata.n_obs}")
print(f"Variables (genes): {adata.n_vars}")

# Check data properties
print(f"\nData sparsity: {(adata.X == 0).sum() / adata.X.size:.2%}")
print(f"Data range: {adata.X.min():.2f} to {adata.X.max():.2f}")

# Store raw data
adata.raw = adata.copy()

print("\n" + "="*50)
print("DIFFERENTIAL EXPRESSION ANALYSIS")
print("="*50)

# Get the conditions from metadata (same as DESeq2 analysis)
conditions = metadata["condition"].value_counts().index[:2].tolist()
print(f"Analyzing contrast between: {conditions}")

condition1, condition2 = conditions[0], conditions[1]
print(f"Condition 1 ({condition1}): {(metadata['condition'] == condition1).sum()} samples")
print(f"Condition 2 ({condition2}): {(metadata['condition'] == condition2).sum()} samples")

# Subset to only include samples with these two conditions
mask = adata.obs['condition'].isin(conditions)
adata_subset = adata[mask].copy()

print(f"Subset data shape: {adata_subset.shape}")

# Method 1: Wilcoxon rank-sum test
print(f"\n--- Wilcoxon Test: {condition1} vs {condition2} ---")
sc.tl.rank_genes_groups(
    adata_subset,
    'condition',
    groups=[condition1],
    reference=condition2,
    method='wilcoxon',
    use_raw=True,
    key_added='wilcoxon_de'
)

# Extract results
wilcox_results = sc.get.rank_genes_groups_df(
    adata_subset,
    group=condition1,
    key='wilcoxon_de',
    pval_cutoff=0.05,
    log2fc_min=0.25
)

print(f"Significant genes (Wilcoxon): {len(wilcox_results)}")
if len(wilcox_results) > 0:
    print("Top 10 upregulated genes:")
    print(wilcox_results.head(10)[['names', 'logfoldchanges', 'pvals_adj']])

# Method 2: t-test
print(f"\n--- T-test: {condition1} vs {condition2} ---")
sc.tl.rank_genes_groups(
    adata_subset,
    'condition',
    groups=[condition1],
    reference=condition2,
    method='t-test',
    use_raw=True,
    key_added='ttest_de'
)

ttest_results = sc.get.rank_genes_groups_df(
    adata_subset,
    group=condition1,
    key='ttest_de',
    pval_cutoff=0.05,
    log2fc_min=0.25
)

print(f"Significant genes (t-test): {len(ttest_results)}")
if len(ttest_results) > 0:
    print("Top 10 upregulated genes:")
    print(ttest_results.head(10)[['names', 'logfoldchanges', 'pvals_adj']])

# Compare methods
if len(wilcox_results) > 0 and len(ttest_results) > 0:
    wilcox_genes = set(wilcox_results['names'])
    ttest_genes = set(ttest_results['names'])
    
    plt.figure(figsize=(8, 6))
    venn2([wilcox_genes, ttest_genes], ('Wilcoxon', 'T-test'))
    plt.title(f'DE Method Comparison: {condition1} vs {condition2}')
    plt.show()
    
    # Overlap statistics
    overlap = wilcox_genes & ttest_genes
    print(f"\nMethod comparison:")
    print(f"Wilcoxon only: {len(wilcox_genes - ttest_genes)}")
    print(f"T-test only: {len(ttest_genes - wilcox_genes)}")
    print(f"Both methods: {len(overlap)}")

# Create volcano plot
print("\n--- Creating Volcano Plot ---")
if len(wilcox_results) > 0:
    # Get all results (not just significant ones)
    all_wilcox_results = sc.get.rank_genes_groups_df(
        adata_subset,
        group=condition1,
        key='wilcoxon_de'
    )
    
    plt.figure(figsize=(10, 6))
    
    # Create volcano plot
    x = all_wilcox_results['logfoldchanges']
    y = -np.log10(all_wilcox_results['pvals_adj'])
    
    # Color points by significance
    colors = ['red' if (abs(fc) > 0.25 and pval < 0.05) else 'gray' 
              for fc, pval in zip(all_wilcox_results['logfoldchanges'], 
                                all_wilcox_results['pvals_adj'])]
    
    plt.scatter(x, y, c=colors, alpha=0.6, s=1)
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 Adjusted P-value')
    plt.title(f'Volcano Plot: {condition1} vs {condition2}')
    
    # Add significance thresholds
    plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    plt.axvline(x=0.25, color='black', linestyle='--', alpha=0.5)
    plt.axvline(x=-0.25, color='black', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.show()

# Save results for comparison with DESeq2
print("\n--- Saving Results ---")
if len(wilcox_results) > 0:
    # Save scanpy results in similar format to DESeq2
    scanpy_results = all_wilcox_results.copy()
    scanpy_results.columns = ['gene', 'log2FoldChange', 'pvalue', 'padj', 'stat']
    scanpy_results = scanpy_results[['gene', 'log2FoldChange', 'stat', 'pvalue', 'padj']]
    scanpy_results.set_index('gene', inplace=True)
    
    # Save to file
    output_path = "./output_files/synthetic_example/"
    import os
    os.makedirs(output_path, exist_ok=True)
    scanpy_results.to_csv(os.path.join(output_path, "scanpy_results.csv"))
    print(f"Scanpy results saved to: {output_path}scanpy_results.csv")

print("\n" + "="*50)
print("COMPARISON WITH DESEQ2 RESULTS")
print("="*50)

# Load DESeq2 results if available
try:
    deseq2_results = pd.read_csv("./output_files/synthetic_example/results.csv", index_col=0)
    print("DESeq2 results loaded successfully!")
    
    # Compare significant genes
    deseq2_sig = deseq2_results[deseq2_results['padj'] < 0.05].index
    scanpy_sig = wilcox_results['names'] if len(wilcox_results) > 0 else []
    
    print(f"DESeq2 significant genes: {len(deseq2_sig)}")
    print(f"Scanpy significant genes: {len(scanpy_sig)}")
    
    if len(deseq2_sig) > 0 and len(scanpy_sig) > 0:
        overlap = set(deseq2_sig) & set(scanpy_sig)
        print(f"Overlapping significant genes: {len(overlap)}")
        
        # Venn diagram comparison
        plt.figure(figsize=(8, 6))
        venn2([set(deseq2_sig), set(scanpy_sig)], ('DESeq2', 'Scanpy'))
        plt.title('Significant Genes: DESeq2 vs Scanpy')
        plt.show()
        
        if len(overlap) > 0:
            print("Overlapping genes:", list(overlap)[:10])  # Show first 10
            
except FileNotFoundError:
    print("DESeq2 results not found. Run DESeq2 analysis first for comparison.")

print("\nAnalysis complete!")
print("You can now compare DESeq2 and scanpy results on the same preprocessed data.")
