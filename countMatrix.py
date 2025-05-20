import os
import sys
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import time

# Configuration
FASTQ_FILE = "data/SRR1552444.fastq"  # Path to your FASTQ file
OUTPUT_DIR = "simple_count_output"
SAMPLE_NAME = "sample1"

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"Processing FASTQ file: {FASTQ_FILE}")

# Check if the FASTQ file exists
if not os.path.exists(FASTQ_FILE):
    print(f"Error: FASTQ file not found at {FASTQ_FILE}")
    sys.exit(1)

# Simple k-mer based approach to count "genes"
# This is a very simplified approach for educational purposes
def process_fastq_to_counts():
    start_time = time.time()
    
    # Define k-mer size
    k = 25  # This would normally be much more sophisticated
    
    # Dictionary to store k-mer counts
    kmer_counts = {}
    
    # Counter for read processing
    read_count = 0
    
    # Open and parse the FASTQ file
    print("Reading FASTQ file and counting k-mers...")
    
    # Determine if the file is gzipped
    is_gzipped = FASTQ_FILE.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    with open_func(FASTQ_FILE, mode) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            read_count += 1
            
            # Print progress every 10,000 reads
            if read_count % 10000 == 0:
                print(f"Processed {read_count} reads...")
            
            # Get the sequence as a string
            seq_str = str(record.seq)
            
            # Generate k-mers and count them
            for i in range(len(seq_str) - k + 1):
                kmer = seq_str[i:i+k]
                kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
            
            # Limit to 100,000 reads for quicker processing
            if read_count >= 100000:
                break
    
    # Convert counts to a dataframe
    print("Converting counts to a dataframe...")
    gene_ids = list(kmer_counts.keys())
    counts = list(kmer_counts.values())
    
    df = pd.DataFrame({
        'gene_id': gene_ids,
        'count': counts
    })
    
    # Sort by count in descending order
    df = df.sort_values('count', ascending=False)
    
    # Keep only the top 20,000 most abundant k-mers (to simulate genes)
    df = df.head(20000)
    
    # Save to CSV
    output_file = os.path.join(OUTPUT_DIR, "gene_counts.csv")
    df.to_csv(output_file, index=False)
    
    print(f"Created count matrix with {len(df)} 'genes'")
    print(f"Saved to {output_file}")
    
    # Create AnnData object
    try:
        import scanpy as sc
        import anndata as ad
        
        # Create AnnData object
        X = np.array(df['count']).reshape(-1, 1)
        var = pd.DataFrame(index=df['gene_id'])
        obs = pd.DataFrame(index=[SAMPLE_NAME])
        
        adata = ad.AnnData(X=X, var=var, obs=obs)
        
        # Add condition for diffxpy
        adata.obs['condition'] = 'condition1'
        
        # Save AnnData object
        h5ad_file = os.path.join(OUTPUT_DIR, "counts.h5ad")
        adata.write_h5ad(h5ad_file)
        print(f"Saved AnnData object to {h5ad_file}")
        
    except ImportError:
        print("Warning: scanpy or anndata not installed. AnnData object not created.")
    
    # Create simple visualization
    plt.figure(figsize=(10, 6))
    plt.hist(np.log10(df['count'] + 1), bins=50)
    plt.xlabel('log10(count + 1)')
    plt.ylabel('Number of k-mers')
    plt.title('Distribution of k-mer Counts')
    plt.savefig(os.path.join(OUTPUT_DIR, "count_distribution.png"))
    plt.close()
    
    # Show top k-mers
    plt.figure(figsize=(12, 8))
    top_df = df.head(20)
    plt.barh(top_df['gene_id'], top_df['count'])
    plt.xlabel('Count')
    plt.ylabel('k-mer')
    plt.title('Top 20 Most Abundant k-mers')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "top_kmers.png"))
    plt.close()
    
    elapsed_time = time.time() - start_time
    print(f"Processing completed in {elapsed_time:.2f} seconds")
    
    return df

# Run the processing
gene_counts = process_fastq_to_counts()

print("\nProcess complete!")
print("Note: This is a simplified approach using k-mers as proxies for genes.")
print("For real research, proper alignment tools like kallisto are recommended.")

input("\nPress Enter to exit...")