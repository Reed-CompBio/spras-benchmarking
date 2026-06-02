import pandas
from pathlib import Path
import collections
from typing import OrderedDict

# from sampling import attempt_sample
import argparse
import networkx

import itertools

import random

panther_directory = Path("..") / ".." / ".."
intermediate_directory = panther_directory / "intermediate"
processed_directory = panther_directory / "processed"

def create_threshold_parser():
    parser = argparse.ArgumentParser(description="Generate thresholded interactomes")
    parser.add_argument("pathway", help="Pathway name (file safe coded)")
    parser.add_argument("--threshold", type=int, required=True, help="Threshold value (0 - 100)")
    # parser.add_argument("--seed", type=int, help="Random seed", default=1234)
    return parser

def count_weights() -> OrderedDict[int, int]:
    """
    Returns an ordered map (lowest to highest weight) from the weight to the number of elements the weight has
    """

    weight_counts = pandas.read_csv(processed_directory / "weight-counts.tsv", sep="\t")
    return collections.OrderedDict(sorted({int(k * 1000): int(v) for k, v in dict(weight_counts.values).items()}.items()))


def find_connected_sources_targets(sources: list[str], targets: list[str], graph: networkx.DiGraph) -> list[tuple[str, str]]:
    connections: list[tuple[str, str]] = []
    for source, target in itertools.product(sources, targets):
        if graph.has_node(source) and graph.has_node(target) and networkx.has_path(graph, source, target):
            connections.append((source, target))
    return connections

def check_target_reachability(graph, sources, targets):
    """
    Check what percentage of targets can be reached from any source
    Returns: (reachability_percentage, reachable_count, total_targets)
    """
    reachable_targets = set()
    for source in sources:
        if graph.has_node(source):
            # Get all nodes reachable FROM this source
            reachable = networkx.descendants(graph, source)
            reachable_targets.update(reachable)
    reachable_count = len(reachable_targets & set(targets))
    total_targets = len(targets)
    reachability = reachable_count / total_targets

    return reachability, reachable_count, total_targets

def sample_interactome(interactome_df: pandas.DataFrame, weight_mapping: OrderedDict[int, int], percentage: float):
    """
    Samples X% of an interactome using its weight_counts dictionary.
    (See `count_weights` for generating `weight_counts`.)

    # Sampling a percentage of the edges from each weight bucket is equivalent to
    # sampling a percentage of the full interactome such that the weight
    # distribution is preserved, since the buckets partition edges by weight.

    percentage is the threshold
    """
    if percentage == 100:
        return interactome_df

    full_list: list[int] = []
    curr_v = 0

    for k, v in weight_mapping.items():
        full_list.extend(map(lambda x: x + curr_v, random.sample(range(1, v), round((percentage/100) * v))))
        curr_v += v
    full_set = set(full_list)

    subset_df = interactome_df.iloc[list(full_set)].reset_index(drop=True)

    return subset_df

def attempt_sample(pathway_df, sampled_interactome, sources, targets, threshold_connectivity=0.8):
    """
    Try to merge pathway with sampled interactome and check if connectivity threshold is met.
    Returns merged_df if successful, None if below threshold.
    """
    median_weight = sampled_interactome["Weight"].median()
    pathway_df["Weight"] = median_weight

    # original = networkx.from_pandas_edgelist(pathway_df, source="Node1", target="Node2", create_using=networkx.DiGraph)
    # original_connections = find_connected_sources_targets(sources, targets, original)

    # merge pathway with sampled interactome
    merged_df = sampled_interactome.merge(
        pathway_df[["Node1", "Node2"]],
        how="left",
        on=["Node1", "Node2"]
    )
    print(len(merged_df))

    # prioritize directed edges over undirected edges, then heightest weight
    merged_df = merged_df.sort_values(
        by=["Node1", "Node2", "Direction", "Weight"],
        ascending=[True, True, False, False]  # D comes before U, highest weight first
    )
    merged_df = merged_df.drop_duplicates(subset=['Node1', 'Node2'], keep='first')

    merged_graph = networkx.from_pandas_edgelist(merged_df, source="Node1", target="Node2", create_using=networkx.DiGraph)
    # merged_connections = find_connected_sources_targets(sources, targets, merged_graph)

    # We ask that at least `percentage` of the sources and targets are connected with one another.
    # connection_percentage = float(len(merged_connections)) / float(len(original_connections)) if len(original_connections) != 0 else 0

    # if threshold_connectivity <= connection_percentage:
    #     print(f"Got {connection_percentage * 100:.1f}% connections above the {threshold_connectivity * 100:.1f}% required percentage threshold.")
    #     return merged_df
    # print(f"Failed {connection_percentage * 100:.1f}% connections below the {threshold_connectivity * 100:.1f}% required percentage threshold.")
    # return None


    # Check reachability
    reachability, reachable_count, total_targets = check_target_reachability(merged_graph, sources, targets)
    print(f"Target reachability: {reachable_count}/{total_targets} ({reachability*100:.1f}%)")

    if reachability >= threshold_connectivity:
        print("Success!")
        return merged_df

    print(f"Below threshold ({threshold_connectivity*100:.0f}%), retrying...")
    return None


def main():
    arg_parser = create_threshold_parser()

    args = arg_parser.parse_args()

    pathway = args.pathway
    # pathway_name = urllib.parse.unquote(pathway)
    # if args.seed is not None:
    #     random.seed(args.seed)

    print(pathway)
    print(args.threshold)

    interactome_df = pandas.read_csv(
        processed_directory / "interactome.tsv",
        header=None,
        sep="\t",
        names=["Node1", "Node2", "Weight", "Direction"],
    )
    pathway_df = pandas.read_csv(
        intermediate_directory / pathway / "pathway.csv",
        sep="\t",
    )

    pathway_input_nodes_df = pandas.read_csv(
        intermediate_directory / pathway / "input_nodes.csv",
        sep="\t",
    )

    weight_mapping = count_weights()

    sources = list(pathway_input_nodes_df[pathway_input_nodes_df["sources"] == True]["NODEID"])
    targets = list(pathway_input_nodes_df[pathway_input_nodes_df["targets"] == True]["NODEID"])

    # if empty sources or targets, then no interactome can be really be used
    if not sources or not targets:
        output_dir = processed_directory / pathway / "interactomes"
        output_dir.mkdir(parents=True, exist_ok=True)

        empty_df = pandas.DataFrame()
        output_path = output_dir / f"interactome_{args.threshold}.txt"
        empty_df.to_csv(output_path, sep="\t", index=False, header=False)

        print(f"No sources or targets found. Wrote empty file to {output_path}")
        return

    # Retry loop
    CONNECTIVITY_THRESHOLD = 0.5  # 80%
    attempt = 1
    merged_df = None
    while True:
        print(f"Attempt {attempt}")
        sampled_interactome = sample_interactome(interactome_df, weight_mapping, args.threshold)

        merged_df = attempt_sample(pathway_df, sampled_interactome, sources, targets, CONNECTIVITY_THRESHOLD)
        if merged_df is not None:
            break
        attempt += 1
        if attempt == 100:
            output_dir = processed_directory / pathway / "interactomes"
            output_dir.mkdir(parents=True, exist_ok=True)

            empty_df = pandas.DataFrame()
            output_path = output_dir / f"interactome_{args.threshold}.txt"
            empty_df.to_csv(output_path, sep="\t", index=False, header=False)

            print(f"Hit 100 attempts. Wrote empty file to {output_path}")
            return

    # Save output
    output_dir = processed_directory / pathway / "interactomes"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"interactome_{args.threshold}.txt"
    merged_df.to_csv(output_path, sep="\t", index=False, header=False)
    print(f"Success on attempt {attempt - 1}")

if __name__ == "__main__":
    main()
