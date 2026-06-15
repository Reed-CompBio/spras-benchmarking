import pandas
from pathlib import Path
import collections
from typing import OrderedDict

import argparse
import networkx

import random

# the directories for the pathways and full interactome
panther_directory = Path("..") / ".." / ".."
intermediate_directory = panther_directory / "intermediate"
processed_directory = panther_directory / "processed"

# directory to put the processed downsampled interactomes
threshold_processed_directory = Path("processed")

# all the checks are done on directed graphs
# all the rest of the data will keep the undirected and directed edges when dealing with combining the pathways and interactomes

def create_threshold_parser():
    parser = argparse.ArgumentParser(description="Generate thresholded interactomes")
    parser.add_argument("pathway", help="Pathway name (file safe coded)")
    parser.add_argument("--threshold", type=int, required=True, help="Threshold value (0 - 100)")
    return parser

def count_weights() -> OrderedDict[int, int]:
    """
    Returns an ordered map (lowest to highest weight) from the weight to the number of elements the weight has
    """

    weight_counts = pandas.read_csv(processed_directory / "weight-counts.tsv", sep="\t")
    return collections.OrderedDict(sorted({int(k * 1000): int(v) for k, v in dict(weight_counts.values).items()}.items()))


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

def to_directed_edgelist(df: pandas.DataFrame) -> pandas.DataFrame:
    directed = df[df["Direction"] == "D"]
    undirected = df[df["Direction"] == "U"]
    flipped = undirected.copy()
    flipped["Node1"], flipped["Node2"] = undirected["Node2"], undirected["Node1"]
    expanded = pandas.concat([directed, undirected, flipped], ignore_index=True)

    expanded = expanded.sort_values("Weight", ascending=False)
    expanded = expanded.drop_duplicates(subset=["Node1", "Node2"], keep="first")
    return expanded.reset_index(drop=True)

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

def combine_pathway_interactome(pathway_df: pandas.DataFrame, sampled_interactome: pandas.DataFrame) -> pandas.DataFrame:
    """
    Concat pathway edges (at median interactome weight) into the sampled interactome.
    Directed edges take priority over undirected; higher weights are kept for duplicates.
    """
    pathway_edges = pathway_df[["Node1", "Node2", "Direction"]].copy()
    pathway_edges["Weight"] = sampled_interactome["Weight"].median()
    pathway_edges = pathway_edges[["Node1", "Node2", "Weight", "Direction"]]

    combined = pandas.concat([sampled_interactome, pathway_edges], ignore_index=True)
    combined = combined.sort_values(
        by=["Node1", "Node2", "Direction", "Weight"],
        ascending=[True, True, True, False],
    )
    return combined.drop_duplicates(subset=["Node1", "Node2"], keep="first").reset_index(drop=True)


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

def check_sample(pathway_df, sampled_interactome, sources, targets, threshold_connectivity=0.8, threshold_overlap=0.25):
    """
    Try to merge pathway with sampled interactome and check if connectivity threshold is met.
    Returns combined df if successful, None if below threshold.
    """

    overlap = pathway_overlap(pathway_df, sampled_interactome)
    if overlap < threshold_overlap:
        print(f"Failed {overlap * 100:.1f}% pathway overlap below the {threshold_overlap * 100:.1f}% threshold.")
        return None, overlap

    print(f"Got {overlap * 100:.1f}% pathway overlap above the {threshold_overlap * 100:.1f}% threshold.")

    # combine pathway and sampled interactome
    combined = combine_pathway_interactome(pathway_df, sampled_interactome)

    # build a directed edge only graph of the combined interactome just for the property checks; U -> two bi-directed D edges
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
        return combined, overlap

    print(f"Below threshold ({threshold_connectivity*100:.0f}%), retrying...")
    return None, overlap

def write_overlap_record(pathway: str, threshold: int, overlap: float):
    """
    Append a single row to interactome-pathway-overlap.tsv.
    Columns: pathway, threshold, overlap
    """
    overlap_path = threshold_processed_directory / pathway / "interactomes" / f"overlap_{threshold}.tsv"
    row = pandas.DataFrame([{"pathway": pathway, "threshold": threshold, "overlap": overlap}])
    row.to_csv(overlap_path, sep="\t", index=False)

def main():
    arg_parser = create_threshold_parser()
    args = arg_parser.parse_args()
    pathway = args.pathway

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
        output_dir = threshold_processed_directory / pathway / "interactomes"
        output_dir.mkdir(parents=True, exist_ok=True)

        empty_df = pandas.DataFrame()
        output_path = output_dir / f"interactome_{args.threshold}.txt"
        empty_df.to_csv(output_path, sep="\t", index=False, header=False)

        write_overlap_record(pathway, args.threshold, 0)

        print(f"No sources or targets found. Wrote empty file to {output_path}")
        return


    # Retry loop
    CONNECTIVITY_THRESHOLD = 1 # 100%
    OVERLAP_THRESHOLD = 0.30 # 30%
    attempt = 1
    MAX_ATTEMPTS = 50
    combined_df = None # the overall best and passed the properties
    best_overlap = 0.0 # this is the best overlap of all the attempts
    best_combined_df = None # this is the best combined interactome of all the attempts

    while True:

        print(f"Attempt {attempt}")
        sampled_interactome = sample_interactome(interactome_df, weight_mapping, args.threshold)

        combined_df, current_overlap = check_sample(pathway_df, sampled_interactome, sources, targets, CONNECTIVITY_THRESHOLD, OVERLAP_THRESHOLD)

        # found a successful combined interactome that passes all the properties
        if combined_df is not None:
            best_overlap = current_overlap
            best_combined_df = combined_df
            break

        # saving the current best incase we don't find an interactome that passes all the properties
        if combined_df is None and current_overlap > best_overlap:
            best_overlap = current_overlap
            best_combined_df = combine_pathway_interactome(pathway_df, sampled_interactome)


        attempt += 1

        if attempt == MAX_ATTEMPTS:
            output_dir = threshold_processed_directory / pathway / "interactomes"
            output_dir.mkdir(parents=True, exist_ok=True)
            output_path = output_dir / f"interactome_{args.threshold}.txt"

            if best_combined_df is not None: # save best attempt
                best_combined_df.to_csv(output_path, sep="\t", index=False, header=False)
                print(f"Hit {MAX_ATTEMPTS} attempts. Wrote best overlap ({best_overlap*100:.1f}%) attempt to {output_path}")
            else:
                pandas.DataFrame().to_csv(output_path, sep="\t", index=False, header=False)
                print(f"Hit {MAX_ATTEMPTS} attempts. Wrote empty file to {output_path}")

            write_overlap_record(pathway, args.threshold, best_overlap)
            return

    # Save successful output
    output_dir = threshold_processed_directory / pathway / "interactomes"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"interactome_{args.threshold}.txt"

    best_combined_df.to_csv(output_path, sep="\t", index=False, header=False)
    write_overlap_record(pathway, args.threshold, best_overlap)
    print(f"Success on attempt {attempt}")

if __name__ == "__main__":
    main()
