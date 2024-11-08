from os import write
from pathlib import Path
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

scores_path = snakemake.input.scores
output_1_path = snakemake.output[0]
output_2_path = snakemake.output[1]

data = []
algorithms = []
pathways = []
valid_pairs = []

for idx, file in enumerate(scores_path):
    file_name = file.split("/")[-1].split(".")[0]
    parts = file_name.split("_")
    if parts[0] not in algorithms:
        algorithms.append(parts[0])
    if parts[1] not in pathways:
        pathways.append(parts[1])
    valid_pairs.append((parts[0], parts[1]))
    with open(file, "r") as f:
        next(f)
        for line in f:
            cols = line.split("\t")
            data.append((parts[0], parts[1], cols[4], cols[5]))

jaccard_edges_indices_list = []
jaccard_nodes_indices_list = []

for algorithm in algorithms:
    current_edge_list = []
    current_node_list = []
    for pathway in pathways:
        appended = False
        for entry in data:
            if entry[0] == algorithm and entry[1] == pathway:
                appended = True
                current_edge_list.append(float(entry[2]))
                current_node_list.append(float(entry[3]))
        if appended == False:
            current_edge_list.append(np.nan)
            current_node_list.append(np.nan)

    jaccard_edges_indices_list.append(current_edge_list)
    jaccard_nodes_indices_list.append(current_node_list)


jaccard_edges_indices = np.array(jaccard_edges_indices_list, dtype=float)
jaccard_nodes_indices = np.array(jaccard_nodes_indices_list, dtype=float)

plt.figure(figsize=(10, 8))
sns.heatmap(
    jaccard_edges_indices,
    annot=True,
    cmap="viridis",
    xticklabels=pathways,
    yticklabels=algorithms,
)

plt.xlabel("Pathways")
plt.ylabel("Algorithms")
plt.title("Jaccard Index Edge Heatmap")
plt.savefig(output_1_path, format="png", dpi=300)

plt.figure(figsize=(10, 8))
sns.heatmap(
    jaccard_nodes_indices,
    annot=True,
    cmap="viridis",
    xticklabels=pathways,
    yticklabels=algorithms,
)

plt.xlabel("Pathways")
plt.ylabel("Algorithms")
plt.title("Jaccard Index Node Heatmap")
plt.savefig(output_2_path, format="png", dpi=300)
