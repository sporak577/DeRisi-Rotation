from Bio import SeqIO

file1 = 'file_1.fasta'
file2 = "file_2.fasta"
output_file = "unique_sequences_by_sequence_only.fasta"

# Load all sequences from file2 into a set
file2_sequences = {str(record.seq) for record in SeqIO.parse(file2, "fasta")}

# Output records from file1 whose sequences are NOT in file2
output_records = []
for record in SeqIO.parse(file1, "fasta"):
    seq = str(record.seq)
    if seq not in file2_sequences:
        output_records.append(record)  # Keep original ID and description

# Write the result
SeqIO.write(output_records, output_file, "fasta")
print(f"Wrote {len(output_records)} unique sequences to {output_file}")
