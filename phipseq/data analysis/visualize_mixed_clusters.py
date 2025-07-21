"""
"""

import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 

# set paths 
csv_file = 'clusters_with_mixed_strains_050825.csv'

#load cluster summary 
df = pd.read_csv(csv_file)

"""
To-DOs 
show distribution arenaviridae vs rest of viral seqs
"""

# -------- PLOT 1: Bar chart of exact cluster sizes --------
cluster_size_counts = df['num_tiles'].value_counts().sort_index()

plt.figure(figsize=(10, 6))
sns.barplot(x=cluster_size_counts.index, y=cluster_size_counts.values)
plt.title("Exact Distribution of Cluster Sizes")
plt.xlabel("Number of Tiles in Cluster")
plt.ylabel("Number of Clusters")
plt.grid(True, axis='y')
plt.tight_layout()
plt.savefig("cluster_size_barplot_exact.png")
plt.close()



# -------- PLOT 3: Top clustered tiles --------
# Expand tile metadata



# -------- Suggested future plots for diversity --------
print("""
Suggested future plots:
- Shannon entropy per tile (computed from AA sequences) vs tile number or position.
- Line plot of entropy across tile numbers to highlight diversity-rich regions.
- Use scatter plot with tile_number (or inferred protein position) on x-axis and entropy on y-axis.
- Overlay protein name or strain color to compare trends across viruses.
- Use a violin or box plot to compare entropy distributions between cluster types.""")
