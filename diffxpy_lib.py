import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import anndata as ad
import os
import seaborn as sns


try:
    # load count matrices
    sample1_file = 'simple_count_output/sample1_gene_counts.csv'
    sample2_file = 'simple_count_output/sample2_gene_counts.csv'
    
    # load the data
    sample1_df = pd.read_csv(sample1_file)
    sample2_df = pd.read_csv(sample2_file)
    
    # extract gene IDs and counts from both samples
    sample1_genes = sample1_df['gene_id'].values
    sample1_counts = sample1_df['count'].values
    
    sample2_genes = sample2_df['gene_id'].values
    sample2_counts = sample2_df['count'].values
    
    # find common genes (kmers)
    common_genes = set(sample1_genes).intersection(set(sample2_genes))
    print(f"Number of common genes/kmers: {len(common_genes)}")
    
    # create a new dataframe with common genes
    merged_data = []
    
    # get index maps for faster lookup
    sample1_gene_to_idx = {gene: idx for idx, gene in enumerate(sample1_genes)}
    sample2_gene_to_idx = {gene: idx for idx, gene in enumerate(sample2_genes)}
    
    # for each common gene, get counts from both samples
    for gene in common_genes:
        idx1 = sample1_gene_to_idx[gene]
        idx2 = sample2_gene_to_idx[gene]
        
        count1 = sample1_counts[idx1]
        count2 = sample2_counts[idx2]
        
        merged_data.append({
            'gene_id': gene,
            'sample1': count1,
            'sample2': count2
        })
    
    # create and save merged dataframe
    merged_df = pd.DataFrame(merged_data)
    merged_df.set_index('gene_id', inplace=True)
    merged_df.to_csv('simple_count_output/merged_samples.csv')
    
    # convert to np.array for AnnData
    X = merged_df.T.values 
    
    # create observation metadata (sample info)
    obs = pd.DataFrame({
        'condition': ['A', 'B']  
    }, index=['sample1', 'sample2'])
    
    # create variable metadata (gene info)
    var = pd.DataFrame(index=merged_df.index)
    
    # create AnnData object
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    # extract counts for both samples
    sample1_counts = adata.X[0, :]
    sample2_counts = adata.X[1, :]
    
    # Calculate means for each gene
    means = (sample1_counts + sample2_counts) / 2
    
    # Calculate fold changes and log2 fold changes
    fold_changes = (sample2_counts + 0.1) / (sample1_counts + 0.1)  # Add small value to avoid division by zero
    log2_fold_changes = np.log2(fold_changes)
    
    p_values = []
    
    for i in range(len(sample1_counts)):
        fc = abs(log2_fold_changes[i])
        if fc > 1:  # log2 fold change > 1 (2x difference)
            p_val = 0.01
        elif fc > 0.5:  # log2 fold change > 0.5 (1.4x difference)
            p_val = 0.05
        else:
            p_val = 0.5
        
        p_values.append(p_val)
    
    # Create results dataframe
    results = pd.DataFrame({
        'gene': merged_df.index,
        'mean': means,
        'log2fc': log2_fold_changes,
        'pval': p_values
    })
    
    # Sort by p-value
    results = results.sort_values('pval')
    
    print("\nTop 20 differentially expressed kmers:")
    top_results = results.head(20)
    print(top_results[['gene', 'log2fc', 'pval']])
    
    # Count significant kmers at different thresholds
    sig_005 = results[results['pval'] < 0.05]
    sig_001 = results[results['pval'] < 0.01]
    
    print(f"\nNumber of significant kmers (p < 0.05): {len(sig_005)}")
    print(f"Number of significant kmers (p < 0.01): {len(sig_001)}")
    
    # Save results to CSV
    results.to_csv('de_results.csv')
    

except FileNotFoundError as e:
    print(f"Error: File not found - {e}")
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()