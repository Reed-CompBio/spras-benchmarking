import argparse
import csv
from pathlib import Path
import sys
import numpy as np
import time
import matplotlib.pyplot as plt
import seaborn as sns
import math
import networkx as nx
from collections import deque, defaultdict
import os


def read_file(path, delimit, type):
    print("Reading", type)
    start_time = time.time()
    with open(path, "r") as file:
        csvreader = csv.reader(file, delimiter=delimit)
        next(csvreader)
        edge_set = set()
        node_set = set()
        edge_confidence_dict = {}
        # for i in range(100):
        #     row = next(csvreader)
        for row in csvreader:
            node_set.add(row[0])
            node_set.add(row[1])
            edge_set.add(frozenset((row[0], row[1])))
            if type == "interactome":
                edge_confidence_dict[frozenset((row[0], row[1]))] = int(row[9])
            else:
                edge_confidence_dict[frozenset((row[0], row[1]))] = int(1000)

    end_time = time.time() - start_time
    print("end time ", end_time)
    return edge_set, node_set, edge_confidence_dict


def get_overlap_stats(human_edge_list, pathway_edge_list):
    count = 0
    for edge in pathway_edge_list:
        if edge not in human_edge_list:
            count += 1
        protein_1 = list(edge)[0]
        protein_2 = list(edge)[1]
    print("Missing edges: ", count)
    return None


def largest_connected_component(edge_set, pathway_edge_list):
    G = nx.Graph()
    G.add_edges_from(map(tuple, edge_set))

    largest_cc = max(nx.connected_components(G), key=len)
    S = G.subgraph(largest_cc)
    print("LCC edges : ", len(S.edges()), " Nodes ", len(S.nodes()))

    lcc_edges = {frozenset(edge) for edge in S.edges()}
    count = sum(1 for edge in pathway_edge_list if edge not in lcc_edges)

    print("Missing", count, "pathway edges out of ", len(pathway_edge_list))


def normalize_array(arr):
    min_val = np.min(arr)
    max_val = np.max(arr)
    return (arr - min_val) / (max_val - min_val)


def main(pathway_path, pathway_name, output_path):

    os.makedirs(output_path, exist_ok=True)
    os.makedirs(f"{output_path}/images", exist_ok=True)
    os.makedirs(f"{output_path}/edge_files", exist_ok=True)
    os.makedirs(f"{output_path}/summary", exist_ok=True)


    human_ppi_path = Path(
        "preprocessing/human-interactome/9606.protein.links.full.v12.0.txt"
    )
    human_edge_set, human_node_set, human_edge_confidence_dict = read_file(
        human_ppi_path, " ", "interactome"
    )

    # pathway = "apoptosis"
    # pathway_path = Path("preprocessing/Apoptosis_signaling/STRING-EDGES-Apoptosis_signaling_.txt")

    pathway_edge_set, pathway_node_set, pathway_edge_confidence_dict = read_file(
        pathway_path, "\t", "pathway"
    )

    # combine both edge lists
    combined_edge_list = human_edge_set.union(pathway_edge_set)
    human_edge_confidence_dict.update(pathway_edge_confidence_dict)

    combined_conf_dict = human_edge_confidence_dict
    print(
        "human interactome edges",
        len(human_edge_set),
        "\npathway interactome edges",
        len(pathway_edge_set),
        "\ntotal interactome edges",
        len(combined_edge_list),
    )

    # stats on missing nodes/edges
    # get_overlap_stats(human_edge_set, pathway_edge_set)

    score_edge_dict = defaultdict(set)

    sorted_combined_conf_dict = {
        key: value
        for key, value in sorted(combined_conf_dict.items(), key=lambda item: item[1])
    }

    for key in sorted_combined_conf_dict:
        score_edge_dict[sorted_combined_conf_dict[key]].add(key)

    edges_removed_count = 0
    pathway_edges_removed_count = 0

    pathway_edges_data = []
    interactome_edges_data = []
    threshold_data = []
    for key in score_edge_dict:
        edges_removed_count += len(score_edge_dict[key])
        pathway_edges_removed_count += len(
            score_edge_dict[key].intersection(pathway_edge_set)
        )

        threshold_data.append(int(key))
        pathway_edges_data.append(len(pathway_edge_set) - pathway_edges_removed_count)
        interactome_edges_data.append(len(combined_edge_list) - edges_removed_count)

        print(
            "threshold:",
            key,
            "\n edges in threshold:",
            len(score_edge_dict),
            "\n pathway edges in threshold:",
            len(score_edge_dict[key].intersection(pathway_edge_set)),
            "\n total edges removed counter:",
            edges_removed_count,
            "\n total pathway edges removed counter:",
            pathway_edges_removed_count,
        )

    results_dict = {
        "alpha_value": [],
        "best_score_threshold": [],
        "interactome_edge_size": [],
        "pathway_edge_size": [],
    }
    results_headers = [
        "alpha_value",
        "best_score_threshold",
        "interactome_edge_size",
        "pathway_edge_size",
    ]
    pathway_arr_norm = normalize_array(np.array(pathway_edges_data))
    interactome_arr_norm = normalize_array(np.array(interactome_edges_data))
    alpha_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    for alpha in alpha_list:
        scores = alpha * pathway_arr_norm - (1 - alpha) * interactome_arr_norm
        best_idx = np.argmax(scores)
        best_threshold = threshold_data[best_idx]

        results_dict["alpha_value"].append(alpha)
        results_dict["best_score_threshold"].append(best_threshold)
        results_dict["interactome_edge_size"].append(interactome_edges_data[best_idx])
        results_dict["pathway_edge_size"].append(pathway_edges_data[best_idx])

        plt.figure(figsize=(14, 5))
        plt.plot(
            threshold_data,
            pathway_arr_norm,
            label="Normalized Pathway Edges",
            marker="o",
            markersize=2,
        )
        plt.plot(
            threshold_data,
            interactome_arr_norm,
            label="Normalized Interactome Edges",
            marker="s",
            markersize=2,
        )
        plt.plot(
            threshold_data,
            scores,
            label="Optimization Score",
            marker="x",
            linestyle="--",
            color="red",
            markersize=2,
        )
        plt.axvline(
            best_threshold,
            color="gray",
            linestyle=":",
            label=f"Optimal Threshold: {best_threshold:.2f}",
        )
        plt.xlabel("Threshold Score")
        plt.ylabel("Normalized Values")
        plt.legend()
        plt.title("Threshold Optimization for Retaining Pathway Edges")
        plt.savefig(Path(f"{output_path}/images/{pathway_name}_{alpha}.pdf"))
        
        edge_file_headers = ["protein1", "protein2", "score"]
        with open(f"{output_path}/edge_files/{pathway_name}_{alpha}.txt", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(edge_file_headers)

            for i in range(best_threshold, 1001):
                if len(score_edge_dict[i]) != 0:
                    edges = score_edge_dict[i]
                    for edge in edges:
                        edge = list(edge)
                        protein_1 = edge[0]
                        protein_2 = edge[1]
                        writer.writerow([protein_1, protein_2, i])
            f.close()


    with open(f"{output_path}/summary/{pathway_name}.csv", "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(results_headers)

        for i in range(len(results_dict["alpha_value"])):
            writer.writerow(
                [
                    results_dict["alpha_value"][i],
                    results_dict["best_score_threshold"][i],
                    results_dict["interactome_edge_size"][i],
                    results_dict["pathway_edge_size"][i],
                ]
            )
        f.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="SPRAS pathway benchmarking network preprocessing"
    )
    parser.add_argument(
        "-p",
        "--pathway_path",
        type=str,
        help="Path to the pathway edge file.",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--pathway_name",
        type=str,
        help="Name of pathway.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        help="Path to the output directory.",
        required=True,
    )

    args = parser.parse_args()

    main(args.pathway_path, args.pathway_name, args.output_dir)
