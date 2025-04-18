#!/bin/bash

# Output file
OUTPUT="combined_proteins_with_accession.faa"
> "$OUTPUT"  # Clear if exists

# Loop through all folders starting with GC*_
for dir in GC*_*; do
    if [ -d "$dir" ] && [ -f "$dir/protein.faa" ]; then
        accession="$dir"
        
        # Modify each header line to include accession ID
        awk -v acc="$accession" '
        /^>/ {print ">" acc "|" substr($0, 2)} 
        !/^>/ {print}
        ' "$dir/protein.faa" >> "$OUTPUT"
    fi
done

echo "Combined protein FASTA written to $OUTPUT"
