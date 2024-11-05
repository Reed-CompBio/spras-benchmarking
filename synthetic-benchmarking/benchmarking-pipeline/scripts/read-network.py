from pathlib import Path
from matplotlib import pyplot as plt
import pandas as pd
import networkx as nx

spras_path = snakemake.input[0]
panther_path = snakemake.input[1]
output_1_path = snakemake.output[0]
output_2_path = snakemake.output[1]


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
print("jaccard_index: ", jaccard_index)

# Calculate Edge Overlap Ratios
overlap_ratio_set1 = (len(intersection) / len(G_set)) * 100
overlap_ratio_set2 = (len(intersection) / len(H_set)) * 100

print("Edge Overlap Ratio (Panther to SPRAS):", overlap_ratio_set1, "%")
print("Edge Overlap Ratio (SPRAS to Panther):", overlap_ratio_set2, "%")


plt.figure(figsize=(10, 8))
pos = nx.spring_layout(G)

nx.draw(
    G,
    pos,
    with_labels=True,
    node_color="lightblue",
    node_size=700,
    font_size=10,
    font_color="black",
    edge_color="gray",
)

plt.title("SPRAS")
plt.savefig(output_1_path, format="png", dpi=300)
plt.close()

plt.figure(figsize=(10, 8))
pos = nx.spring_layout(H)

nx.draw(
    H,
    pos,
    with_labels=True,
    node_color="lightblue",
    node_size=700,
    font_size=10,
    font_color="black",
    edge_color="gray",
)

plt.title("PANTHER")
plt.savefig(output_2_path, format="png", dpi=300)
plt.close()
