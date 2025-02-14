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
