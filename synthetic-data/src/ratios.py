import pandas as pd
import os

directory = "spras-compatible-pathway-data/"

folders = [f for f in os.listdir(directory)]

ratio = pd.DataFrame(columns=['pathway', '#_edges', '#_targets', '#_sources', 'ratio_targets', 'ratio_sources'])
for f in folders:
    edges_file = f"{directory}{f}/{f}_gs_edges.txt"
    node_prizes_file = f"{directory}{f}/{f}_node_prizes.txt"

    edges_df = pd.read_csv(edges_file, sep = "\t", header=None)
    node_prizes_df = pd.read_csv(node_prizes_file, sep = "\t", header=0)


    num_edges = len(edges_df)
    num_sources = node_prizes_df["sources"].sum()
    num_targets = node_prizes_df["targets"].sum()

    ratio_sources = round(num_sources / num_edges, 3)
    ratio_targets = round(num_targets / num_edges, 3)

    new_row = {
        'pathway': f,
        '#_edges': num_edges,
        '#_targets': num_targets,
        '#_sources': num_sources,
        'ratio_targets': ratio_targets,
        'ratio_sources': ratio_sources,
    }   

    ratio = pd.concat([ratio, pd.DataFrame([new_row])], ignore_index=True)  

ratio.to_csv("spras-compatible-pathway-data/data_ratios.txt", sep="\t", index = False)