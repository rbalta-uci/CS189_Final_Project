import os
import gzip
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

# input and output files
FASTQ_FILES = [
    {"file": "data/SRR1552444.fastq", "sample": "sample1"},
    {"file": "data/SRR1552445.fastq", "sample": "sample2"} 
]
OUTPUT_DIR = "count_output_virgin"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# process fastq file and return k-mer counts
def fastq_to_counts(fastq, sample, max_reads=100000, k=25):
    
    kmer_counts = {}
    read_count = 0

    is_gzipped = fastq.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    with open_func(fastq, mode) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            read_count += 1
            
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
    
    # top 20000 k-mers
    df = df.head(20000)

    # save count files
    output_file = os.path.join(OUTPUT_DIR, f"{sample}_gene_counts.csv")
    df.to_csv(output_file, index=False)
    
    return df

# process FASTQ files
all_samples_data = {}
merged_kmers = set()

for sample_info in FASTQ_FILES:
    fastq = sample_info["file"]
    sample = sample_info["sample"]
    
    df = fastq_to_counts(fastq, sample)
    if df is not None:
        all_samples_data[sample] = df