"""
Authors: ChatGPT and Sophie-Christine Porak
"""
import os
import pandas as pd
from glob import glob

# === Set your directory path ===
data_dir = 'your_directory'
file_paths = glob(os.path.join(data_dir, "*.csv"))  # change to *.tsv if needed

# === Initialize empty list for dataframes ===
dfs = []

for file in file_paths:
    sample_id = os.path.splitext(os.path.basename(file))[0]  # get sample ID from filename

    try:
        # Read only the first two columns
        df = pd.read_csv(file, usecols=[0, 1], header=None, names=["peptide", sample_id])

        # Coerce non-numeric counts to NaN and drop those rows
        df[sample_id] = pd.to_numeric(df[sample_id], errors='coerce')
        df = df.dropna(subset=[sample_id])

        dfs.append(df.set_index("peptide"))
    
    except Exception as e:
        print(f"Skipping {file}: {e}")

# === Merge all sample columns on Peptide ID ===
count_matrix = pd.concat(dfs, axis=1).fillna(0).astype(int)

# === Save to file (optional) ===
count_matrix.to_csv("peptide_count_matrix.csv")

# === View result ===
print(count_matrix.head())
