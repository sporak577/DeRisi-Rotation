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

date = "051225.02"

# Function to make pie chart with 3% cutoff
"""
series is a pandas seroies of category labels,
label is a string used to plot title and file naming,
outname is a filename stem for the output image,
weights is optional series of numeric values,
cutoff is the proportion threshold"""
def plot_pie(series, label, outname, weights=None, cutoff=0.02):
    if weights is not None:
        # Group by category and for each group it sums the weights, here the tile counts. gives total count per category. 
        counts = series.groupby(series).apply(lambda x: weights.loc[x.index].sum())
    else:
        counts = series.value_counts()
    #if no weights are provided is simply counts the number of times each category appears in the series
    total = counts.sum()
    #calculates the total number of tiles (or total weigh)
    filtered = counts[counts / total >= cutoff]  # e.g., 3% threshold, filters out all categories that contribute less than a specific cutoff
    other_sum = total - filtered.sum() #calculates how many counts were excluded, will be added as a single other category
    #if any values were excluded, add a new entry called Other (<2%) to the filtered results
    if other_sum > 0:
        filtered[f"Other (<{int(cutoff * 100)}%)"] = other_sum

    plt.figure(figsize=(14, 14))
    plt.pie(
        filtered,
        labels=filtered.index,
        autopct='%1.1f%%',
        startangle=140,
        colors=colors[:len(filtered)],
        textprops={'fontsize': 15}
    )
    plt.title(f"Tile count representation of viral {label}", size=15)
    plt.tight_layout()
    plt.savefig(f"/Users/sophieporak/Desktop/{outname}_{date}.png", dpi=300)
    plt.close()


# Apply to different taxonomic levels   
plot_pie(df['tax_rank_7'], "Families", "virus_family_distribution_pie", weights=df['count'])
plot_pie(df['tax_rank_8'], "Genera", "virus_genus_distribution_pie", weights=df['count'])
plot_pie(df['organism'], "Strains", "virus_strain_distribution_pie", weights=df['count'])
plot_pie(df['protein_name'], "Proteins", "virus_protein_distribution_pie", weights=df['count'])


# -------- Filtering only for Mammarenavirus genus --------

#filter for only mammarenavirus genus 
df.columns = df.columns.str.strip()
mammarenavirus_df = df[df['tax_rank_8'] == 'Mammarenavirus']

print(mammarenavirus_df['protein_name'].value_counts())

#normalizing synonymous protein names
protein_rename_map = {
    "L protein": "L protein",
    "polymerase": "L protein",
    "RNA-directed RNA polymerase": "L protein",
    "RNA dependent-RNA polymerase": "L protein",
    "RNA-dependent RNA polymerase": "L protein",
    "L polymerase": "L protein",
    "L": "L protein",
    "RNA-directed RNA polymerase L": "L protein",
    "RNA polymerase": "L protein",
    "large RNA-dependent RNA polymerase": "L protein",
    "large RNA-dependent RNA polymerase protein": "L protein",
    "RNA dependent RNA polymerase": "L protein",
    "polymerase RDRP": "L protein",

    "hypothetical protein": "Z protein",
    "Z protein": "Z protein",
    "Z-protein": "Z protein",
    "RING finger protein Z": "Z protein",
    "ring finger protein": "Z protein",
    "small RING finger protein": "Z protein",
    "zinc binding protein": "Z protein",
    "zinc finger protein": "Z protein",
    "zinc finger-like protein": "Z protein",
    "RING finger Z protein": "Z protein",
    "zinc-binding protein": "Z protein",
    "multifunctional matrix-like protein": "Z protein",
    "matrix protein": "Z protein",

    "glycoprotein precursor": "GPC",
    "glycoprotein": "GPC",
    "preglycoprotein polyprotein GP complex": "GPC",
    "glycoprotein precursor complex": "GPC",
    "Glycoprotein precursor": "GPC",
    "glycoprotein G1+G2 precursor": "GPC",
    "envelope glycoprotein": "GPC",
    "GPC": "GPC",
    "GP": "GPC",

    "nucleocapsid protein": "NP",
    "nucleoprotein": "NP",
    "NP": "NP",
    "nucleocapsid": "NP",
    

}

plot_pie(
    mammarenavirus_df['protein_name'].replace(protein_rename_map),
    "Proteins within Mammarenavirus genus",
    "mammarenavirus_protein_tile_distribution_pie", 
    weights=mammarenavirus_df['count']
)


# -------- Filtering only for Mammarenavirus lassaense --------

df.columns = df.columns.str.strip()
mammarenavirus_df = df[df['organism'] == 'Mammarenavirus lassaense']

print(mammarenavirus_df['protein_name'].value_counts())


plot_pie(
    mammarenavirus_df['protein_name'].replace(protein_rename_map),
    "Proteins within Mammarenavirus lassaense",
    "mammarenavirus_lassaense_protein_tile_distribution_pie",
    weights=mammarenavirus_df['count']
)

