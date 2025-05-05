from Bio import SeqIO
import re

file1 = "file1.fasta"
file2 = "file2.fasta"
output_file = "unique_sequences_annotated.fasta"

def extract_root_id(record):
    """Gets the root ID like XOB76279.1 from the FASTA header."""
    return record.id.split("_")[0]

def extract_protein_and_tile(description):
    """
    Extracts protein name and tile number from full description.
    Example: 'nucleoprotein|...| tile 7 of 22' → ('nucleoprotein', '7')
    """
    parts = description.split("|")
    protein = parts[0].split()[-1] if parts else "unknown"

    # Look for 'tile N of M'
    tile_match = re.search(r'tile\s+(\d+)\s+of\s+\d+', description)
    tile = tile_match.group(1) if tile_match else "?"

    return protein, tile

# Load all root IDs from file2
file2_ids = {extract_root_id(record) for record in SeqIO.parse(file2, "fasta")}

# Compare and extract annotated records from file1
output_records = []
for record in SeqIO.parse(file1, "fasta"):
    root_id = extract_root_id(record)
    if root_id not in file2_ids:
        protein, tile = extract_protein_and_tile(record.description)
        new_header = f"{root_id} | {protein} | tile {tile}"
        record.id = new_header
        record.description = ""
        output_records.append(record)

# Write output
SeqIO.write(output_records, output_file, "fasta")
print(f"Wrote {len(output_records)} unique annotated records to {output_file}")
