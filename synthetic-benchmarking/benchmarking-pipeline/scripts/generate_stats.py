from pathlib import Path
from matplotlib import pyplot as plt
import pandas as pd
import networkx as nx

spras_path = snakemake.input[0]
panther_path = snakemake.input[1]
output_path = snakemake.output[0]

spras_df = pd.read_csv(spras_path, sep="\t")
panther_df = pd.read_csv(panther_path, sep="\t")

G = nx.from_pandas_edgelist(spras_df, source="Node1", target="Node2")
H = nx.from_pandas_edgelist(panther_df, source="Node1", target="Node2")

G_set = set(G.edges())
H_set = set(H.edges())
union = G_set | H_set
intersection = G_set & H_set
jaccard_index = len(intersection) / len(union)

print("union: ", union)
print("intersection: ", intersection)
print("Jaccard Index: ", jaccard_index)

# Calculate Edge Overlap Ratios
overlap_ratio_set1 = (len(intersection) / len(G_set)) * 100
overlap_ratio_set2 = (len(intersection) / len(H_set)) * 100

f = open(output_path, "w+")
f.write(f"Jaccard Index: {jaccard_index}\n")
f.write(f"Edge Overlap Ratio (Panther to SPRAS): {overlap_ratio_set1}%\n")
f.write(f"Edge Overlap Ratio (SPRAS to Panther):: {overlap_ratio_set2}%\n")
