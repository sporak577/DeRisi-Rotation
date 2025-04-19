"""
this script extracts information from a full genebank file 
and outputs a fasta file with the header information you define
to be extracted. 
"""

from Bio import SeqIO

# === Settings ===
input_file = "your_file.gb"    # Replace with your .gb or .gbff file
output_file = "proteins_with_metadata.fasta"

with open(output_file, "w") as out_f:
    for record in SeqIO.parse(input_file, "genbank"):
        accession = record.id
        # Get metadata from source feature
        source = next((f for f in record.features if f.type == "source"), None)
        geo = source.qualifiers.get("geo_loc_name", ["?"])[0] if source else "?"
        host = source.qualifiers.get("host", ["?"])[0] if source else "?"
        date = source.qualifiers.get("collection_date", ["?"])[0] if source else "?"

        # Loop through CDS features that have translation
        for feature in record.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                product = feature.qualifiers.get("product", ["?"])[0]
                aa_seq = feature.qualifiers["translation"][0]

                header = f">{accession} | {product} | {geo} | {host} | {date}"
                out_f.write(header + "\n")
                # Wrap sequence every 70 characters
                for i in range(0, len(aa_seq), 70):
                    out_f.write(aa_seq[i:i+70] + "\n")

print(f"Done! Protein FASTA saved to: {output_file}")
