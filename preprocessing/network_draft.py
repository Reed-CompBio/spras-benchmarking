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


def read_file(path, delimit, type, prior_knowledge_score):
    print("Reading", type)
    start_time = time.time()
    with open(path, "r") as file:
        csvreader = csv.reader(file, delimiter=delimit)
        next(csvreader)
        edge_set = set()
        node_set = set()
        edge_confidence_dict = {}
        # for i in range(3000):
        #     row = next(csvreader)
        for row in csvreader:
            node_set.add(row[0])
            node_set.add(row[1])
            edge_set.add(frozenset((row[0], row[1])))
            if type == "interactome":
                edge_confidence_dict[frozenset((row[0], row[1]))] = int(row[9])
            else:
                edge_confidence_dict[frozenset((row[0], row[1]))] = (
                    prior_knowledge_score
                )

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

def combine_pathway_files(pathway_dir, pathway_file_list, output_dir):

    edge_set = set()
    for pathway_file in pathway_file_list:
        print(f"{pathway_dir}/{pathway_file}")
        with open(f"{pathway_dir}/{pathway_file}", "r") as f:
            csvreader = csv.reader(f, delimiter="\t")
            next(csvreader)
            for row in csvreader:
                edge_set.add(frozenset((row[0], row[1])))
        f.close()

    with open(f"{output_dir}/edge_files/combined_pathway_edge_file.txt", "w") as f:
        edge_file_headers = ["protein1", "protein2"]
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(edge_file_headers)

        for edge in edge_set:
            edge = list(edge)
            protein_1 = edge[0]
            protein_2 = edge[1]
            writer.writerow([protein_1, protein_2])
        f.close()
    
    return None


def main(pathway_dir, output_path):

    os.makedirs(output_path, exist_ok=True)
    os.makedirs(f"{output_path}/edge_files", exist_ok=True)
    os.makedirs(f"{output_path}/summary", exist_ok=True)

    combine_pathway_files(pathway_dir, os.listdir(pathway_dir), output_path)


    human_ppi_path = Path(
        "preprocessing/human-interactome/9606.protein.links.full.v12.0.txt"
    )
    human_edge_set, human_node_set, human_edge_confidence_dict = read_file(
        human_ppi_path, " ", "interactome", None
    )

    sorted_human_edge_score_dict = {
        key: value
        for key, value in sorted(
            human_edge_confidence_dict.items(), key=lambda item: item[1]
        )
    }

    results = {
        "score_threshold": [],
        "edges_in_threshold_human_interacome": [],
        "edges_in_pathway": [],
        "added_pathway_edges": [],
        "edges_in_combined_filtered_interacome": [],
    }

    for threshold in range(1, 1001, 100):

        score_edge_dict = defaultdict(set)

        threhold_filtered_human_edge_count = 0
        score_dist = []
        for key in sorted_human_edge_score_dict:
            current_score = human_edge_confidence_dict[key]
            if current_score >= int(threshold):
                score_edge_dict[current_score].add(key)
                score_dist.append(current_score)
                threhold_filtered_human_edge_count += 1

        percentile = 50
        score_dist_arr = np.array(score_dist)
        prior_knowledge_score = int(np.percentile(score_dist_arr, percentile))
        print(prior_knowledge_score)

        pathway_edge_set, pathway_node_set, pathway_edge_confidence_dict = read_file(
            f"{output_path}/edge_files/combined_pathway_edge_file.txt", "\t", "pathway", prior_knowledge_score
        )

        added_pathway_edges_count = 0
        for key in pathway_edge_confidence_dict:
            current_score = pathway_edge_confidence_dict[key]
            if key not in human_edge_confidence_dict:
                score_edge_dict[current_score].add(key)
                added_pathway_edges_count += 1
            elif human_edge_confidence_dict[key] < int(threshold):
                score_edge_dict[current_score].add(key)
                added_pathway_edges_count += 1

        combined_edge_list = set()

        for score in score_edge_dict:
            combined_edge_list.update(score_edge_dict[score])

        print("number of edges in original human interactome", len(human_edge_set))
        print(
            "number of edges in threhold filtered human interactome",
            threhold_filtered_human_edge_count,
        )
        print("number of pathway edges", len(pathway_edge_set))
        print("number of total added pathway edges", added_pathway_edges_count)
        print("number of total combined filtered interacome", len(combined_edge_list))

        edge_file_headers = ["protein1", "protein2", "score"]
        with open(f"{output_path}/edge_files/combined_{threshold}.txt", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(edge_file_headers)

            for score in score_edge_dict:
                for edge in score_edge_dict[score]:
                    edge = list(edge)
                    protein_1 = edge[0]
                    protein_2 = edge[1]
                    writer.writerow([protein_1, protein_2, score])

            f.close()

        results["score_threshold"].append(threshold)
        results["edges_in_threshold_human_interacome"].append(
            threhold_filtered_human_edge_count
        )
        results["edges_in_pathway"].append(len(pathway_edge_set))
        results["added_pathway_edges"].append(added_pathway_edges_count)
        results["edges_in_combined_filtered_interacome"].append(len(combined_edge_list))

    results_headers = [
        "score_threshold",
        "edges_in_threshold_human_interacome",
        "edges_in_pathway",
        "added_pathway_edges",
        "edges_in_combined_filtered_interacome",
    ]
    with open(f"{output_path}/summary/combined.csv", "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(results_headers)

        for i in range(len(results["score_threshold"])):
            writer.writerow(
                [
                    results["score_threshold"][i],
                    results["edges_in_threshold_human_interacome"][i],
                    results["edges_in_pathway"][i],
                    results["added_pathway_edges"][i],
                    results["edges_in_combined_filtered_interacome"][i],
                ]
            )
        f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="SPRAS pathway benchmarking network preprocessing"
    )
    parser.add_argument(
        "-p",
        "--pathway_dir",
        type=str,
        help="Path to the pathway edge file.",
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

    main(args.pathway_dir, args.output_dir)
