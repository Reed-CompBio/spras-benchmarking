import pandas
from pathlib import Path
import collections
from typing import OrderedDict

import argparse
import networkx

import random

panther_directory = Path("..") / ".." / ".."
intermediate_directory = panther_directory / "intermediate"
processed_directory = panther_directory / "processed"

# all the checks are done on directed graphs
# all the rest of the data will keep the undirected and directed edges

def create_threshold_parser():
    parser = argparse.ArgumentParser(description="Generate thresholded interactomes")
    parser.add_argument("pathway", help="Pathway name (file safe coded)")
    parser.add_argument("--threshold", type=int, required=True, help="Threshold value (0 - 100)")
    # parser.add_argument("--seed", type=int, help="Random seed", default=1234) # TODO add a seed
    return parser

def count_weights() -> OrderedDict[int, int]:
    """
    Returns an ordered map (lowest to highest weight) from the weight to the number of elements the weight has
    """

    weight_counts = pandas.read_csv(processed_directory / "weight-counts.tsv", sep="\t")
    return collections.OrderedDict(sorted({int(k * 1000): int(v) for k, v in dict(weight_counts.values).items()}.items()))


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
    reachability = reachable_count / total_targets if total_targets else 0.0

    return reachability, reachable_count, total_targets

def check_source_reachability(graph, sources, targets):
    """
    Fraction of sources that can reach at least one target.
    Returns: (reachability, reachable_count, total_sources)
    """
    target_set = set(targets)
    reachable_sources = 0
    for source in sources:
        if graph.has_node(source) and (networkx.descendants(graph, source) & target_set):
            reachable_sources += 1
    total_sources = len(sources)
    reachability = reachable_sources / total_sources if total_sources else 0.0
    return reachability, reachable_sources, total_sources


def to_directed_pairs(df):
    """
    Set of directed (Node1, Node2) edges.
    'D' rows contribute one directed edge. 'U' rows contribute both directions.
    A frame without a Direction column is treated as fully undirected.
    """
    pairs = set(zip(df["Node1"], df["Node2"]))
    if "Direction" in df.columns:
        und = df[df["Direction"] == "U"]
        pairs |= set(zip(und["Node2"], und["Node1"]))
    else:
        pairs |= set(zip(df["Node2"], df["Node1"]))
    return pairs


def pathway_overlap(pathway_df, sampled_interactome):
    """
    Fraction of pathway edges present in the sampled interactome, using the
    directed model the rest of the pipeline uses. A directed pathway edge
    matches a same-direction background edge or an undirected background edge.
    """
    pathway_pairs = to_directed_pairs(pathway_df)
    if not pathway_pairs:
        return 0.0
    sample_pairs = to_directed_pairs(sampled_interactome)
    return len(pathway_pairs & sample_pairs) / len(pathway_pairs)


def sample_interactome(interactome_df: pandas.DataFrame, weight_mapping: OrderedDict[int, int], percentage: float):
    """
    Samples X% of an interactome using its weight_counts dictionary.
    (See `count_weights` for generating `weight_counts`.)

    # Sampling a percentage of the edges from each weight bucket is equivalent to
    # sampling a percentage of the full interactome such that the weight
    # distribution is preserved, since the buckets partition edges by weight.

    percentage is the threshold of how much we take from the full interactome
    """
    if percentage == 100:
        return interactome_df

    full_list: list[int] = []
    curr_v = 0

    for k, v in weight_mapping.items():
        sample_size = min(round((percentage / 100) * v), v)
        full_list.extend(map(lambda x: x + curr_v, random.sample(range(v), sample_size)))
        curr_v += v
    full_set = set(full_list)

    subset_df = interactome_df.iloc[list(full_set)].reset_index(drop=True)

    return subset_df

def to_directed_edgelist(df: pandas.DataFrame) -> pandas.DataFrame:
    directed = df[df["Direction"] == "D"]
    undirected = df[df["Direction"] == "U"]
    flipped = undirected.copy()
    flipped["Node1"], flipped["Node2"] = undirected["Node2"], undirected["Node1"]
    expanded = pandas.concat([directed, undirected, flipped], ignore_index=True)

    expanded = expanded.sort_values("Weight", ascending=False)
    expanded = expanded.drop_duplicates(subset=["Node1", "Node2"], keep="first")
    return expanded.reset_index(drop=True)

def check_sample(pathway_df, sampled_interactome, sources, targets, threshold_connectivity=0.8, threshold_overlap=0.25):
    """
    Try to merge pathway with sampled interactome and check if connectivity threshold is met.
    Returns combined df if successful, None if below threshold.
    """

    overlap = pathway_overlap(pathway_df, sampled_interactome)
    if overlap < threshold_overlap:
        print(f"Failed {overlap * 100:.1f}% pathway overlap below the {threshold_overlap * 100:.1f}% threshold.")
        return None

    print(f"Got {overlap * 100:.1f}% pathway overlap above the {threshold_overlap * 100:.1f}% threshold.")

    # set pathway edges to the median edge weights in the current sampled interactome
    pathway_edges = pathway_df[["Node1", "Node2", "Direction"]].copy()
    pathway_edges["Weight"] = sampled_interactome["Weight"].median()
    pathway_edges = pathway_edges[["Node1", "Node2", "Weight", "Direction"]]

    # concat pathway with sampled interactome
    # the full pathway is added to the sampled interactome
    combined = pandas.concat([sampled_interactome, pathway_edges], ignore_index=True)
    combined = combined.sort_values(
        by=["Node1", "Node2", "Direction", "Weight"],
        ascending=[True, True, True, False],  # want to keep directed edges first (D < U) and then pick highest edge weight for any more duplicate edges
    )
    combined = combined.drop_duplicates(subset=["Node1", "Node2"], keep="first")

    # build a directed edge only graph of the combined interactome just for the checks; U -> two bi-directed D edges
    check_graph = networkx.from_pandas_edgelist(
        to_directed_edgelist(combined),
        source="Node1",
        target="Node2",
        create_using=networkx.DiGraph,
    )

    # Check reachability
    # asking what fraction of targets are reachable from at least one source and what fraction of sources are reachable from at least one target
    reachability_targets, reachable_count_targets, total_targets = check_target_reachability(check_graph, sources, targets)
    reachability_sources, reachable_count_sources, total_sources = check_source_reachability(check_graph, sources, targets)
    print(f"Target reachability: {reachable_count_targets}/{total_targets} ({reachability_targets*100:.1f}%)")
    print(f"Source reachability: {reachable_count_sources}/{total_sources} ({reachability_sources*100:.1f}%)")

    if reachability_targets >= threshold_connectivity and reachability_sources >= threshold_connectivity:
        print("Success!")
        return combined

    print(f"Below threshold ({threshold_connectivity*100:.0f}%), retrying...")
    return None


def main():
    arg_parser = create_threshold_parser()

    args = arg_parser.parse_args()

    pathway = args.pathway
    # if args.seed is not None:
    #     random.seed(args.seed)

    print(f"Pathway: {pathway}")
    print(f"Thresholding by: {args.threshold}")

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
        output_dir = Path("processed") / pathway / "interactomes"
        output_dir.mkdir(parents=True, exist_ok=True)

        empty_df = pandas.DataFrame()
        output_path = output_dir / f"interactome_{args.threshold}.txt"
        empty_df.to_csv(output_path, sep="\t", index=False, header=False)

        print(f"No sources or targets found. Wrote empty file to {output_path}")
        return

    # Retry loop
    CONNECTIVITY_THRESHOLD = 1
    OVERLAP_THRESHOLD = 0.30 # 50% might be too high, maybe we do like 45% or 40%? or go back to 25%
    attempt = 1
    combined_df = None

    while True:
        print(f"Attempt {attempt}")
        sampled_interactome = sample_interactome(interactome_df, weight_mapping, args.threshold)

        combined_df = check_sample(pathway_df, sampled_interactome, sources, targets, CONNECTIVITY_THRESHOLD, OVERLAP_THRESHOLD)
        if combined_df is not None:
            break

        attempt += 1

        if attempt == 100:
            # Save unsuccessful output
            output_dir = Path("processed") / pathway / "interactomes"
            output_dir.mkdir(parents=True, exist_ok=True)

            empty_df = pandas.DataFrame()
            output_path = output_dir / f"interactome_{args.threshold}.txt"
            empty_df.to_csv(output_path, sep="\t", index=False, header=False)

            print(f"Hit 100 attempts. Wrote empty file to {output_path}")
            return

    # Save successful output
    output_dir = Path("processed") / pathway / "interactomes"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"interactome_{args.threshold}.txt"
    combined_df.to_csv(output_path, sep="\t", index=False, header=False)
    print(f"Success on attempt {attempt}")

if __name__ == "__main__":
    main()
