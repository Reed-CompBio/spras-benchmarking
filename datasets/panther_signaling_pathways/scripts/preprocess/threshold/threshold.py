import pandas
from pathlib import Path
import collections
from typing import OrderedDict

# from sampling import attempt_sample
import argparse
import networkx

import random

# TODO: think about directionality
# - i think treat everything as a directed edge when doing any of the checks
# - but do keep things as undirected and directed as everything else

panther_directory = Path("..") / ".." / ".."
intermediate_directory = panther_directory / "intermediate"
processed_directory = panther_directory / "processed"


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

def pathway_overlap(pathway_df, sampled_interactome):
    """
    Fraction of pathway edges present in the sampled interactome, ignoring direction.
    Edges are compared as unordered pairs.
    """
    pathway_pairs = {
        frozenset((a, b)) for a, b in zip(pathway_df["Node1"], pathway_df["Node2"])
    }
    if not pathway_pairs:
        return 0.0

    sample_pairs = {
        frozenset((a, b))
        for a, b in zip(sampled_interactome["Node1"], sampled_interactome["Node2"])
    }
    return len(pathway_pairs & sample_pairs) / len(pathway_pairs)

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

def to_directed_edgelist(df: pandas.DataFrame) -> pandas.DataFrame:
    directed = df[df["Direction"] == "D"]
    undirected = df[df["Direction"] == "U"]
    flipped = undirected.copy()
    flipped["Node1"], flipped["Node2"] = undirected["Node2"], undirected["Node1"]
    return pandas.concat([directed, undirected, flipped], ignore_index=True)

def attempt_sample(pathway_df, sampled_interactome, sources, targets, threshold_connectivity=0.8, threshold_overlap=0.25):
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
    pathway_edges = pathway_df[["Node1", "Node2"]].copy()
    pathway_edges["Weight"] = sampled_interactome["Weight"].median()

     # concat pathway with sampled interactome
    combined = pandas.concat([sampled_interactome, pathway_edges], ignore_index=True)
    combined = combined.sort_values(
        by=["Node1", "Node2", "Direction", "Weight"],
        ascending=[True, True, True, False],  # 'D' < 'U', so ascending keeps directed first, highest weight first as well
    )
    combined = combined.drop_duplicates(subset=["Node1", "Node2"], keep="first")

    combined_graph = networkx.from_pandas_edgelist(combined, source="Node1", target="Node2", create_using=networkx.DiGraph)

    # Check reachability
    # asking what fraction of targets are reachable from at least one source and what fraction of sources are reachable from at least one target
    reachability_targets, reachable_count_targets, total_targets = check_target_reachability(combined_graph, sources, targets)
    reachability_sources, reachable_count_sources, total_sources = check_source_reachability(combined_graph, sources, targets)
    print(f"Target reachability: {reachable_count_targets}/{total_targets} ({reachability_targets*100:.1f}%)")
    print(f"Source reachability: {reachable_count_sources}/{total_sources} ({reachability_sources*100:.1f}%)")

    # TODO: see if I can change the threshold_connectivity to 80% instead of 50%
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
    CONNECTIVITY_THRESHOLD = 0.8
    OVERLAP_THRESHOLD = 0.25
    attempt = 1
    combined_df = None

    while True:
        print(f"Attempt {attempt}")
        sampled_interactome = sample_interactome(interactome_df, weight_mapping, args.threshold)

        combined_df = attempt_sample(pathway_df, sampled_interactome, sources, targets, CONNECTIVITY_THRESHOLD, OVERLAP_THRESHOLD)
        if combined_df is not None:
            break

        attempt += 1

        if attempt == 15:
            # Save unsuccessful output
            output_dir = Path("processed") / pathway / "interactomes"
            output_dir.mkdir(parents=True, exist_ok=True)

            empty_df = pandas.DataFrame()
            output_path = output_dir / f"interactome_{args.threshold}.txt"
            empty_df.to_csv(output_path, sep="\t", index=False, header=False)

            print(f"Hit 15 attempts. Wrote empty file to {output_path}")
            return

    # Save successful output
    output_dir = Path("processed") / pathway / "interactomes"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"interactome_{args.threshold}.txt"
    combined_df.to_csv(output_path, sep="\t", index=False, header=False)
    print(f"Success on attempt {attempt}")

if __name__ == "__main__":
    main()
