"""
This script extracts protein sequences and metadata from a GenPept file and writes a FASTA file
with headers that include accession, product, geographic location, host, collection date, segment, and strain.
The output is sorted by accession number.
"""

from Bio import SeqIO
from operator import itemgetter

print("Biopython works!")

# === Settings ===
input_file = "/Users/sophieporak/Library/CloudStorage/Box-Box/DeRisi/Arenavirus/20250424/arenaviridae_042425.gp"
output_file = "arenaviridae_042425_gb.fasta"

# === Collect all entries first ===
entries = []

for record in SeqIO.parse(input_file, "genbank"): #GenPept uses the same parser
    accession = record.id #gets the accession number from the record
    product = record.description #gets the description line for the record
    aa_seq = str(record.seq) #gets the actual protein sequence from the ORIGIN section

    # Get metadata from source feature
    source = next((f for f in record.features if f.type == "source"), None)
    geo = source.qualifiers.get("geo_loc_name", ["?"])[0] if source else "?"
    host = source.qualifiers.get("host", ["?"])[0] if source else "?"
    date = source.qualifiers.get("collection_date", ["?"])[0] if source else "?"
    segment = source.qualifiers.get("segment", ["?"])[0] if source else "?"
    strain = source.qualifiers.get("strain", ["?"])[0] if source else "?"

   
    entries.append({
            "accession": accession,
            "product": product,
            "geo": geo,
            "host": host,
            "date": date,
            "segment": segment,
            "strain": strain,
            "sequence": aa_seq
    })

# === Sort by accession number ===
entries.sort(key=itemgetter("accession"))

# === Write to output ===
with open(output_file, "w") as out_f:
    for entry in entries:
        header = (
            f">{entry['accession']} | {entry['product']} | {entry['geo']} | "
            f"{entry['host']} | {entry['date']} | segment={entry['segment']} | strain={entry['strain']}"
        )
        out_f.write(header + "\n")
        for i in range(0, len(entry["sequence"]), 70):
            out_f.write(entry["sequence"][i:i+70] + "\n")

print(f"Done! Protein FASTA saved to: {output_file}")
