from Bio import SeqIO
import matplotlib.pyplot as plt
import diffxpy.api as de
import numpy as np

fastq_file = "data/SRR1552444.fastq" 

read_lengths = []
quality_scores = []

max_reads = 100000
count = 0

for record in SeqIO.parse(fastq_file, "fastq"):
    if count >= max_reads:
        break
    
    read_lengths.append(len(record.seq))
    
    if isinstance(record.letter_annotations["phred_quality"][0], int):
        quality_scores.append(np.mean(record.letter_annotations["phred_quality"]))
    else:
        quality_scores.append(np.mean([ord(q) - 33 for q in record.letter_annotations["phred_quality"]]))
    
    count += 1

plt.figure(figsize=(10, 6))
plt.hist(read_lengths, bins=50)
plt.xlabel('Read Length')
plt.ylabel('Count')
plt.title('Read Length Distribution')
plt.savefig('read_length_distribution.png')
plt.show()

plt.figure(figsize=(10, 6))
plt.hist(quality_scores, bins=50)
plt.xlabel('Average Quality Score')
plt.ylabel('Count')
plt.title('Quality Score Distribution')
plt.savefig('quality_score_distribution.png')
plt.show()

print(f"Total reads processed: {count}")
print(f"Average read length: {np.mean(read_lengths):.2f}")
print(f"Average quality score: {np.mean(quality_scores):.2f}")
import diffxpy.api as de

print("hello")