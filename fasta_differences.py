from Bio import SeqIO
import re

file1 = '/Users/sophieporak/Library/CloudStorage/Box-Box/DeRisi/0.94_final_library_protein_records_tiling_out/0.94_tiling_out_aa_preprocess_tiles.fasta'
file2 = "/Users/sophieporak/Library/CloudStorage/Box-Box/DeRisi/cd-hit_0.94_final_library_peptide_records_tiling_out/0.94_tiling_out_aa_preprocess_tiles.fasta"
output_file = "unique_sequences_by_sequence_only.fasta"

def extract_protein_and_tile(description):
    """
    Extracts protein name and tile number from full description.
    Example: 'nucleoprotein|...| tile 7 of 22' → ('nucleoprotein', '7')
    """
    parts = description.split("|")
    protein = parts[0].split()[-1] if parts else "unknown"
    tile_match = re.search(r'tile\s+(\d+)\s+of\s+\d+', description)
    tile = tile_match.group(1) if tile_match else "?"
    return protein, tile

# Load all sequences from file2 into a set (sequence-level comparison)
file2_sequences = {str(record.seq) for record in SeqIO.parse(file2, "fasta")}

# Output records from file1 whose sequences are NOT in file2
output_records = []
for record in SeqIO.parse(file1, "fasta"):
    seq = str(record.seq)
    if seq not in file2_sequences:
        protein, tile = extract_protein_and_tile(record.description)
        record.id = f"{protein} | tile {tile}"
        record.description = ""
        output_records.append(record)

# Write the result
SeqIO.write(output_records, output_file, "fasta")
print(f"Wrote {len(output_records)} unique sequences to {output_file}")
