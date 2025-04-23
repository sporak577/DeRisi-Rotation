import matplotlib.pyplot as plt 
import pandas as pd 
import numpy as np 

thresholds = [0.95, 0.96, 0.97, 0.98, 0.99, 1.0]
unique_tiles = [16275, 17515, 19720, 21851, 23842, 25137]

# convert thresholds to strings so they appear clearly on x-axis
threshold_labels = [str(t) for t in thresholds]

plt.figure(figsize=(8,5))
bars = plt.bar(threshold_labels, unique_tiles)

for bar in bars: 
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, height + 200, f'{height}',
             ha = 'center', va = 'bottom', fontsize = 9)


# Optional grid — helps visually compare heights
plt.grid(axis='y', linestyle='--', alpha=0.7)

plt.xlabel("CD-HIT Identity Threshold")
plt.ylabel("Number of Unique Arenavirus Library Tiles")
plt.show()