from pathlib import Path
from matplotlib import pyplot as plt
import pandas as pd
import networkx as nx

spras_path = snakemake.input[0]
panther_path = snakemake.input[1]
output_path = snakemake.output[0]

spras_df = pd.read_csv(spras_path, sep="\t")
panther_df = pd.read_csv(panther_path, sep="\t")

spras_df = spras_df[["Node1", "Node2"]]
panther_df = panther_df[["Node1", "Node2"]]

all_edges = pd.concat([spras_df, panther_df]).drop_duplicates()
edge_list = set(zip(all_edges["Node1"], all_edges["Node2"]))
ground_truth = {edge: 1 for edge in zip(panther_df["Node1"], panther_df["Node2"])}

print(all_edges)
print(edge_list)
print(ground_truth)

y_true = []
y_scores = []

f = open(output_path, "w+")
f.write("Node1\tNode2\ty_true\ty_score\n")
for edge in edge_list:
    y_true.append(ground_truth.get(edge, 0))  # 1 if edge in file 1, else 0
    y_scores.append(1 if edge in zip(spras_df["Node1"], spras_df["Node2"]) else 0)
    f.write(f"{edge[0]}\t{edge[1]}\t{y_true[-1]}\t{y_scores[-1]}\n")
