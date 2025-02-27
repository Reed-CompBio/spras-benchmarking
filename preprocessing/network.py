import csv
from pathlib import Path
import numpy as np
import time
import matplotlib.pyplot as plt
import seaborn as sns
import math
import networkx as nx
from collections import deque


def read_file(path, delimit, type):
    print("Reading", type)
    start_time = time.time()
    with open(path, "r") as file:
        csvreader = csv.reader(file, delimiter=delimit)
        next(csvreader)
        edge_set = set()
        node_set = set()
        edge_confidence_dict = {}
        # for i in range(10):
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


def main():
    print("Hello")
    human_ppi_path = Path("human-interactome/9606.protein.links.full.v12.0.txt")
    human_edge_list, human_node_set, human_edge_confidence_dict = read_file(
        human_ppi_path, " ", "interactome"
    )

    # pathways = ["hedgehog"]
    pathways = ["apoptosis", "cadherin", "hedgehog", "notch", "wnt"]
    # pathways = ["cadherin", "notch", "wnt"]

    # delimiter = [" ", " ", "\t"]

    pathway_paths = [
        Path("Apoptosis_signaling/STRING-EDGES-Apoptosis_signaling_.txt"),
        Path("Cadherin_signaling/STRING-EDGES-Cadherin_signaling_p.txt"),
        Path("Hedgehog_signaling/STRING-EDGES-Hedgehog_signaling_p.txt"),
        Path("Notch_signaling/STRING-EDGES-Notch_signaling_path.txt"),
        Path("Wnt_signaling/STRING-EDGES-Wnt_signaling_pathwa.txt"),
    ]

    for pathway, pathway_path in zip(pathways, pathway_paths):

        pathway_edge_list, pathway_node_set, pathway_edge_confidence_dict = read_file(
            pathway_path, "\t", "pathway"
        )

        # combine both edge lists
        combined_edge_list = human_edge_list.union(pathway_edge_list)
        combined_conf_dict = pathway_edge_confidence_dict | human_edge_confidence_dict
        print(
            len(human_edge_list),
            len(pathway_edge_list),
            len(combined_edge_list),
            len(combined_conf_dict),
        )

        # stats on missing nodes/edges
        get_overlap_stats(human_edge_list, pathway_edge_list)

        # thresholding
        max_conf_score = 0
        score_dist = []
        min_score = 0
        for key in combined_conf_dict:
            if int(combined_conf_dict[key]) >= min_score:
                max_conf_score = max(max_conf_score, int(combined_conf_dict[key]))
                score_dist.append(int(combined_conf_dict[key]))
        score_dist = sorted(score_dist)
        unique_score_list = sorted(list(set(score_dist)))
        # fig = plt.figure(figsize=(14, 6))
        # plt.hist(score_dist)
        # plt.show()

        # sns.kdeplot(score_dist)
        # plt.show()

        sorted_combined_conf_dict = {
            key: value
            for key, value in sorted(
                combined_conf_dict.items(), key=lambda item: item[1]
            )
        }

        filtered_edge_list = deque(sorted_combined_conf_dict.keys())
        i = 0
        pathway_edges_data = []
        interactome_edges_data = []
        threshold_data = []
        edges_removed = 0
        top_edge = filtered_edge_list.popleft()
        removed_edges_set = set()
        removed_edges_set.add(top_edge)
        for threshold in unique_score_list:
            print("Threshold ", threshold)
            count = 0
            # print(top_edge,  combined_conf_dict[top_edge], threshold, type(combined_conf_dict[top_edge]), type(threshold))
            # print(len(filtered_edge_list) != 0)
            # print(combined_conf_dict[top_edge] == threshold)
            while (
                len(filtered_edge_list) != 0
                and combined_conf_dict[top_edge] == threshold
            ):
                # print("INSIDEEEE")
                # print(filtered_edge_list[0], combined_conf_dict[top_edge])
                top_edge = filtered_edge_list.popleft()
                removed_edges_set.add(top_edge)
                count += 1
                # print(count)
            print("Removed ", count, "edges", "total edges ", len(combined_conf_dict))
            edges_removed += count
            # check largest connected component
            # largest_connected_component(filtered_edge_list, pathway_edge_list)
            print(
                "pathway edges removed : ",
                len(pathway_edge_list.difference(set(removed_edges_set))),
                "/",
                len(pathway_edge_list),
            )
            pathway_edges_data.append(
                len(pathway_edge_list)
                - (
                    len(pathway_edge_list)
                    - len(pathway_edge_list.difference(set(removed_edges_set)))
                )
            )
            threshold_data.append(threshold)
            interactome_edges_data.append(len(combined_conf_dict) - edges_removed)
            print(
                "pathway_edges_data",
                len(pathway_edge_list)
                - len(pathway_edge_list.difference(set(filtered_edge_list))),
            )
            print("threshold_data", threshold)
            print("interactome_edges_data", len(combined_conf_dict) - edges_removed)
            print()

        # print("Threshold", threshold_data)
        # print("interactome edges", interactome_edges_data)
        # print("pathway edges", pathway_edges_data)

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
            # alpha = 0.5
            scores = alpha * pathway_arr_norm - (1 - alpha) * interactome_arr_norm
            best_idx = np.argmax(scores)
            best_threshold = threshold_data[best_idx]

            results_dict["alpha_value"].append(alpha)
            results_dict["best_score_threshold"].append(best_threshold)
            results_dict["interactome_edge_size"].append(
                interactome_edges_data[best_idx]
            )
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
            plt.savefig(Path(f"output/images/{pathway}_{alpha}.pdf"))
            # plt.show()

        with open(f"output/{pathway}.csv", "w") as f:
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


if __name__ == "__main__":
    main()
