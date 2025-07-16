"""
Authors: ChatGPT & Sophie-Christine Porak 
"""

import csv

input_fasta = '0.96_tiling_out_nt_tiles_cp_cleaned_no_linkers.fasta'
output_csv = "parsed_tiles_with_description_and_virus.csv"

records = []

with open(input_fasta) as f:
    seq_id = ref_seq_id = sequence = protein_description = virus_name = None

    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if seq_id and sequence:
                records.append([seq_id, ref_seq_id, protein_description, virus_name, sequence])
            parts = line[1:].split()
            seq_id = parts[0]
            ref_seq_id = parts[2]

            # Rebuild the description string
            description = " ".join(parts[3:])
            # Split by "|" to extract the fields
            fields = description.split("|")
            protein_description = fields[0].strip() if len(fields) > 0 else ""
            virus_name = fields[1].strip() if len(fields) > 1 else ""
            sequence = ""
        else:
            sequence += line

    # Append the last record
    if seq_id and sequence:
        records.append([seq_id, ref_seq_id, protein_description, virus_name, sequence])

# Write to CSV
with open(output_csv, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow(["seq_id", "ref_seq_id", "protein_description", "virus_name", "sequence"])
    writer.writerows(records)

print(f"Saved to {output_csv}")
