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
overlap_ratio_set1 = (len(intersection) / len(G_set)) * 100
overlap_ratio_set2 = (len(intersection) / len(H_set)) * 100

columns = [
    "spras_nodes",
    "panther_nodes",
    "spras_edges",
    "panther_edges",
    "jaccard_index",
    "panther_spras_edge_overlap",
    "spras_panther_edge_overlap",
]

rows = [
    str(len(G.nodes())),
    str(len(H.nodes())),
    str(len(G.edges())),
    str(len(H.edges())),
    str(jaccard_index),
    str(overlap_ratio_set1),
    str(overlap_ratio_set2),
]

f = open(output_path, "w+")
column = "\t".join(columns) + "\t"
row = "\t".join(rows) + "\t"
f.write(f"{column}\n")
f.write(f"{row}\n")