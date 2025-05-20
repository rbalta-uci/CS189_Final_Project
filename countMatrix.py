import os
import sys
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import time

# Define both input files and output directories
FASTQ_FILES = [
    {"file": "data/SRR1552444.fastq", "sample": "sample1"},
    {"file": "data/SRR1552445.fastq", "sample": "sample2"}  # Added second FASTQ file
]
OUTPUT_DIR = "simple_count_output"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Function to process a FASTQ file and return k-mer counts
def process_fastq_to_counts(fastq_file, sample_name, max_reads=100000, k=25):
    start_time = time.time()
    
    print(f"Processing FASTQ file: {fastq_file} as {sample_name}")
    
    if not os.path.exists(fastq_file):
        print(f"Error: FASTQ file not found at {fastq_file}")
        return None
    
    kmer_counts = {}
    read_count = 0

    is_gzipped = fastq_file.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    with open_func(fastq_file, mode) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            read_count += 1
            
            if read_count % 10000 == 0:
                print(f"Processed {read_count} reads...")
            
            seq_str = str(record.seq)
            
            for i in range(len(seq_str) - k + 1):
                kmer = seq_str[i:i+k]
                kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
            
            if read_count >= max_reads:
                break

    gene_ids = list(kmer_counts.keys())
    counts = list(kmer_counts.values())
    
    df = pd.DataFrame({
        'gene_id': gene_ids,
        'count': counts
    })

    df = df.sort_values('count', ascending=False)
    
    # Take top 20000 k-mers
    df = df.head(20000)

    # Save individual count file
    output_file = os.path.join(OUTPUT_DIR, f"{sample_name}_gene_counts.csv")
    df.to_csv(output_file, index=False)
    print(f"Saved counts to {output_file}")
    
    elapsed_time = time.time() - start_time
    print(f"Processed {read_count} reads in {elapsed_time:.2f} seconds")
    
    return df

# Process both FASTQ files
all_samples_data = {}
merged_kmers = set()

for sample_info in FASTQ_FILES:
    fastq_file = sample_info["file"]
    sample_name = sample_info["sample"]
    
    df = process_fastq_to_counts(fastq_file, sample_name)
    if df is not None:
        all_samples_data[sample_name] = df
        merged_kmers.update(df['gene_id'])

# Create a merged count matrix with all k-mers
print(f"Creating merged count matrix with {len(merged_kmers)} unique k-mers")

# Initialize the merged dataframe with the k-mer list
merged_df = pd.DataFrame(index=list(merged_kmers))

# Add each sample's counts as a column
for sample_name, df in all_samples_data.items():
    # Create a series with the sample's counts
    sample_counts = pd.Series(df['count'].values, index=df['gene_id'])
    
    # Add this as a column to the merged dataframe, filling missing values with 0
    merged_df[sample_name] = sample_counts
    
# Fill NaN values with 0
merged_df = merged_df.fillna(0)

# Save the merged matrix
merged_output = os.path.join(OUTPUT_DIR, "merged_count_matrix.csv")
merged_df.to_csv(merged_output)
print(f"Saved merged count matrix to {merged_output}")

# Create an AnnData object from the merged matrix
try:
    import anndata as ad
    
    # Convert to AnnData
    X = merged_df.T.values  # samples x features
    var = pd.DataFrame(index=merged_df.index)
    obs = pd.DataFrame(index=merged_df.columns)
    
    obs['condition'] = ['A', 'B']  
    
    adata = ad.AnnData(X=X, var=var, obs=obs)
    
    h5ad_file = os.path.join(OUTPUT_DIR, "merged_counts.h5ad")
    adata.write_h5ad(h5ad_file)
    print(f"Saved merged AnnData object to {h5ad_file}")
    
except ImportError:
    print("Warning: scanpy or anndata not installed. AnnData object not created.")

print("Analysis complete!")