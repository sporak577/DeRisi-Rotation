"""
differences based on FASTA header - so here the protein IDs
"""

from Bio import SeqIO

file1 = "file1.fasta"
file2 = "file2.fasta"
output_file = "differences_by_id.fasta"

def get_header_id(record):
    """
    Extracts the first part of the FASTA header before the first space or pipe.
    E.g., >XOB76282.1_21_1 -> XOB76282.1
    """
    header = record.id  # takes the part after '>'
    return header.split('_')[0].split('|')[0]  # strip tile info or anything after _ or |

# Load file2 IDs to compare against
file2_ids = {get_header_id(record) for record in SeqIO.parse(file2, "fasta")}

# Compare file1 and output entries not in file2
diff_records = [
    record for record in SeqIO.parse(file1, "fasta")
    if get_header_id(record) not in file2_ids
]

# Write the difference records to output
SeqIO.write(diff_records, output_file, "fasta")
print(f"Wrote {len(diff_records)} records with unique headers to {output_file}")
