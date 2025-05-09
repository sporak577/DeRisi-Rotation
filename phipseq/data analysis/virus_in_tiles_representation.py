"""
Input: 
a CSV file of the metadata of the phage-display tiles that we sent off for ordering. 

this file has following columns: 
accession, count, protein_name, organism, isolate, host, geo_loc_name, 
collection_date, segment, strain, gb_accession, taxnonmy_lineage

Output: 
Pie chart for family, genus and species distribution. 
"""

import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.cm as cm

colors = cm.get_cmap('tab20').colors  # or 'Set3', 'Paired', etc.

df = pd.read_csv('/Users/sophieporak/Library/CloudStorage/Box-Box/DeRisi/data processing/lib_data_AK_050725/0.96_tiling_out_nt_tiles_cp_cleaned_metadata_curated.csv')

date = "050925.0.01"

# Function to make pie chart with 2% cutoff
def plot_pie(series, label, outname):
    counts = series.value_counts() #index is the category and value is the count
    total = counts.sum() #total number of items (sum of all counts)
    filtered = counts[counts / total >= 0.01] #filters counts to only keep labels that represent omore or equal 2% of total 
    other_sum = total - filtered.sum()
    
    if other_sum > 0: #sums how many items were excluded
        filtered["Other (<1%)"] = other_sum

    plt.figure(figsize=(14, 14))
    plt.pie(filtered, labels=filtered.index, autopct='%1.1f%%', startangle=140, colors=colors[:len(filtered)], textprops={'fontsize': 15})
    plt.title(f"Tile count representation of viral {label}", size = 15)
    plt.tight_layout()
    plt.savefig(f"/Users/sophieporak/Desktop/{outname}_{date}.png", dpi = 300)
    plt.close()

# Apply to different taxonomic levels
plot_pie(df['tax_rank_7'], "Families", "virus_family_distribution_pie")
plot_pie(df['tax_rank_8'], "Genera", "virus_genus_distribution_pie")
plot_pie(df['organism'], "Strains", "virus_strain_distribution_pie")
plot_pie(df['protein_name'], "Proteins", "virus_protein_distribution_pie")