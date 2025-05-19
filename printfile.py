def parse_and_print_fastq(fastq_file):
    """
    Parse and print a FASTQ file with its components labeled.
    
    Args:
        fastq_file (str): Path to the FASTQ file
    """
    with open(fastq_file, 'r') as file:
        line_count = 0
        read_count = 1
        
        for line in file:
            line = line.strip()
            if line_count % 4 == 0:
                print(f"Header: {line}")
            elif line_count % 4 == 1:
                print(f"Sequence: {line}")
            elif line_count % 4 == 2:
                print(f"Separator: {line}")
            elif line_count % 4 == 3:
                print(f"Quality: {line}")
                print("-" * 50)  
                read_count += 1  
            
            line_count += 1

# Example usage
fastq_path = "your_file.fastq"
parse_and_print_fastq("data\SRR1552444.fastq")