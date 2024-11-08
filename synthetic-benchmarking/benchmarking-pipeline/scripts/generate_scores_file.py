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

spras_df = spras_df[["Node1", "Node2"]]
panther_df = panther_df[["Node1", "Node2"]]

all_edges = pd.concat([spras_df, panther_df]).drop_duplicates()
edge_list = set(zip(all_edges["Node1"], all_edges["Node2"]))
ground_truth_edge = {edge: 1 for edge in zip(panther_df["Node1"], panther_df["Node2"])}

y_true_edge = []
y_scores_edge = []


f = open(output_1_path, "w+")
f.write("Node1\tNode2\ty_true\ty_score\n")
for edge in edge_list:
    y_true_edge.append(ground_truth_edge.get(edge, 0))  # 1 if edge in file 1, else 0
    y_scores_edge.append(1 if edge in zip(spras_df["Node1"], spras_df["Node2"]) else 0)
    f.write(f"{edge[0]}\t{edge[1]}\t{y_true_edge[-1]}\t{y_scores_edge[-1]}\n")


all_nodes = pd.concat([spras_df, panther_df]).drop_duplicates()
node_list = pd.unique(all_nodes[["Node1", "Node2"]].values.ravel())
spras_node_list = pd.unique(spras_df[["Node1", "Node2"]].values.ravel())
panther_node_list = pd.unique(panther_df[["Node1", "Node2"]].values.ravel())

print(spras_path)
print(f"spras nodes {spras_node_list}")
print(f"spras nodes {panther_node_list}")
print(f"all nodes {node_list}")

ground_truth_node = {node: 1 for node in panther_node_list}
print(ground_truth_node)

y_true_node = []
y_scores_node = []

g = open(output_2_path, "w+")
g.write("Node1\ty_true\ty_score\n")
for node in node_list:
    y_true_node.append(ground_truth_node.get(node, 0))  # 1 if node in file 1, else 0
    y_scores_node.append(1 if node in spras_node_list else 0)
    g.write(f"{node}\t{y_true_node[-1]}\t{y_scores_node[-1]}\n")

