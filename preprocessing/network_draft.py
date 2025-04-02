import argparse
import csv
from pathlib import Path
import re
import sys
import numpy as np
import time
import matplotlib.pyplot as plt
import seaborn as sns
import math
import networkx as nx
from collections import deque, defaultdict
import os
import requests
import time
from requests.adapters import HTTPAdapter, Retry
from urllib.parse import urlparse, parse_qs, urlencode
import zlib
import json
from xml.etree import ElementTree


POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"

retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def read_file(path, delimit, type, prior_knowledge_score):
    print("Reading", type)
    start_time = time.time()
    with open(path, "r") as file:
        csvreader = csv.reader(file, delimiter=delimit)
        next(csvreader)
        edge_set = set()
        node_set = set()
        edge_confidence_dict = {}
        for i in range(3000):
            row = next(csvreader)
        # for row in csvreader:
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


def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    check_response(request)
    return request.json()["jobId"]


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] in ("NEW", "RUNNING"):
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def get_id_mapping_results_search(url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results


def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")


def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


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

    output_results = {
        "score_threshold": [],
        "edges_in_threshold_human_interacome": [],
        "edges_in_pathway": [],
        "added_pathway_edges": [],
        "edges_in_combined_filtered_interacome": [],
    }

    thresholds = [1, 100, 200, 300, 400, 500, 600, 700, 800, 900]

    for threshold in thresholds:

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
            f"{output_path}/edge_files/combined_pathway_edge_file.txt",
            "\t",
            "pathway",
            prior_knowledge_score,
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

        combined_edge_set = set()
        combined_node_set = set()

        for score in score_edge_dict:
            combined_edge_set.update(score_edge_dict[score])
            for edge in score_edge_dict[score]:
                edge = list(edge)
                protein1 = edge[0]
                protein2 = edge[1]
                combined_node_set.add(protein1)
                combined_node_set.add(protein2)

        job_id = submit_id_mapping(
            from_db="STRING", to_db="UniProtKB", ids=list(combined_node_set)
        )

        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link)

            mapped_entries = {}
            counter = 0
            for i in range(len(results["results"])):
                if (
                    results["results"][i]["from"] != ""
                    and results["results"][i]["to"]["primaryAccession"] != ""
                ):
                    mapped_entries[results["results"][i]["from"]] = results["results"][
                        i
                    ]["to"]["primaryAccession"]
                    counter += 1
                else:
                    "NOT MAPPED"
            print("mapped counter ", counter)

            original_ids = combined_node_set
            mapped_ids = set(mapped_entries.keys())
            failed_ids = original_ids - mapped_ids

            print("successful IDs mapped: ", len(mapped_ids))
            print("failed ids mapped: ", len(failed_ids))
            print(list(failed_ids))
            print(
                "pathway edges not mapped ",
                len(failed_ids.intersection(pathway_node_set)),
                "/ ",
                len(failed_ids),
            )

            # remove nodes and edges that does not have a uniprot mapping
            combined_node_set = combined_node_set - failed_ids
            temp_edge_set = set()

            for edge in combined_edge_set:
                edge = list(edge)
                protein_1 = edge[0]
                protein_2 = edge[1]
                if protein_1 not in failed_ids and protein_2 not in failed_ids:
                    temp_edge_set.add(frozenset((protein_1, protein_2)))
                
            combined_edge_set = temp_edge_set

            print("number of edges in original human interactome", len(human_edge_set))
            print(
                "number of edges in threhold filtered human interactome",
                threhold_filtered_human_edge_count,
            )
            print("number of pathway edges", len(pathway_edge_set))
            print("number of total added pathway edges", added_pathway_edges_count)
            print(
                "number of total combined filtered interacome", len(combined_edge_set)
            )

            edge_file_headers = ["protein1", "protein2", "score"]
            with open(f"{output_path}/edge_files/combined_{threshold}.txt", "w") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerow(edge_file_headers)

                for score in score_edge_dict:
                    for edge in score_edge_dict[score]:
                        if edge in combined_edge_set:
                            edge = list(edge)
                            protein_1 = edge[0]
                            protein_2 = edge[1]
                            writer.writerow([mapped_entries[protein_1], mapped_entries[protein_2], score])

                f.close()

            output_results["score_threshold"].append(threshold)
            output_results["edges_in_threshold_human_interacome"].append(
                threhold_filtered_human_edge_count
            )
            output_results["edges_in_pathway"].append(len(pathway_edge_set))
            output_results["added_pathway_edges"].append(added_pathway_edges_count)
            output_results["edges_in_combined_filtered_interacome"].append(
                len(combined_edge_set)
            )

        output_results_headers = [
            "score_threshold",
            "edges_in_threshold_human_interacome",
            "edges_in_pathway",
            "added_pathway_edges",
            "edges_in_combined_filtered_interacome",
        ]
        with open(f"{output_path}/summary/combined.csv", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(output_results_headers)

            for i in range(len(output_results["score_threshold"])):
                writer.writerow(
                    [
                        output_results["score_threshold"][i],
                        output_results["edges_in_threshold_human_interacome"][i],
                        output_results["edges_in_pathway"][i],
                        output_results["added_pathway_edges"][i],
                        output_results["edges_in_combined_filtered_interacome"][i],
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
