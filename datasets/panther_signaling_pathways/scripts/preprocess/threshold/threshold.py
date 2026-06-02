import pandas
from pathlib import Path
import collections
from typing import OrderedDict, NamedTuple

import urllib.parse
from sampling import attempt_sample
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

# def count_weights() -> OrderedDict[int, int]:
#     """
#     Returns an ordered map (lowest to highest weight) from the weight to the number of elements the weight has
#     """

#     weight_counts = pandas.read_csv(processed_directory / "weight-counts.tsv", sep="\t")
#     return collections.OrderedDict(sorted({int(k * 1000): int(v) for k, v in dict(weight_counts.values).items()}.items()))

def count_weights() -> OrderedDict[int, int]:
    weight_counts = pandas.read_csv(processed_directory / "weight-counts.tsv", sep="\t")
    print("weight_counts dataframe:")
    print(weight_counts)
    print("\nweight_counts.values:")
    print(weight_counts.values)
    
    result = collections.OrderedDict(sorted({int(k * 1000): int(v) for k, v in dict(weight_counts.values).items()}.items()))
    print(f"\nweight_mapping has {len(result)} entries")
    return result

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

def attempt_sample(pathway_df, sampled_interactome, sources, targets, original_graph, threshold_connectivity=0.8):
    """
    Try to merge pathway with sampled interactome and check if connectivity threshold is met.
    Returns merged_df if successful, None if below threshold.
    """
    # merge pathway with sampled interactome
    merged_df = pathway_df.merge(
        sampled_interactome[["Node1", "Node2"]], 
        how="inner", 
        on=["Node1", "Node2"]
    )
    # priotize directed edges over undirected edges, then heightest weight
    merged_df = merged_df.sort_values(
        by=["Node1", "Node2", "Direction", "Weight"],
        ascending=[True, True, False, False]  # D comes before U, highest weight first
    )
    merged_df = merged_df.drop_duplicates(subset=['Node1', 'Node2'], keep='first')

    print(f"Merged df shape: {merged_df.shape}")
    
    # Build graph
    graph = networkx.from_pandas_edgelist(
        merged_df, 
        source="Node1", 
        target="Node2", 
        create_using=networkx.DiGraph
    )

    print(f"Graph nodes: {graph.number_of_nodes()}, edges: {graph.number_of_edges()}")
    
    # Check reachability
    reachability, reachable_count, total_targets = check_target_reachability(graph, sources, targets)
    print(f"Target reachability: {reachable_count}/{total_targets} ({reachability*100:.1f}%)")
    
    if reachability >= threshold_connectivity:
        print(f"Success!")
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
    # print(args.seed)
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

    # 1) sample the ineractomes using the weight counts in a stratified manner for a given threshold
  
    weight_mapping = count_weights()
    # sampled_interactome = sample_interactome(interactome_df, weight_mapping, args.threshold)
    
    # 2) take the pathway and merge it into the sampled interactome

    # median_weight = sampled_interactome["Weight"].median()
    # pathway_df["Weight"] = median_weight
   
    # 3) merge the edges of the pathway with what is in the sampled interactome
    # merged_df = pathway_df.merge(
    #     sampled_interactome[["Node1", "Node2"]], 
    #     how="inner", 
    #     on=["Node1", "Node2"]
    # )

    # # priotize directed edges over undirected edges, then heightest weight
    # merged_df = merged_df.sort_values(
    #     by=["Node1", "Node2", "Direction", "Weight"],
    #     ascending=[True, True, False, False]  # D comes before U, highest weight first
    # )
    # merged_df = merged_df.drop_duplicates(subset=['Node1', 'Node2'], keep='first')

    sources = list(pathway_input_nodes_df[pathway_input_nodes_df["sources"] == True]["NODEID"])
    targets = list(pathway_input_nodes_df[pathway_input_nodes_df["targets"] == True]["NODEID"])

    # if empty sources or targets, then no interactome
    if not sources or not targets:
        output_dir = processed_directory / pathway / "interactomes"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        empty_df = pandas.DataFrame()
        output_path = output_dir / f"interactome_{args.threshold}.txt"
        empty_df.to_csv(output_path, sep="\t", index=False, header=False)
        
        print(f"No sources or targets found. Wrote empty file to {output_path}")
        return 

    # Build original graph for reference
    # original_graph = networkx.from_pandas_edgelist(
    #     pathway_df, 
    #     source="Node1", 
    #     target="Node2", 
    #     create_using=networkx.DiGraph
    # )


    # 4) If properties aren't passing, restart the sampling
    # while True:
        
    #     # Build graph from merged edges
    #     graph = networkx.from_pandas_edgelist(
    #         merged_df, 
    #         source="Node1", 
    #         target="Node2", 
    #         create_using=networkx.DiGraph
    #     )
    #     sources_in_graph = [s for s in sources if graph.has_node(s)]
    #     print(f"Sources in graph: {len(sources_in_graph)}/{len(sources)}")
    #     targetes_in_graph = [t for t in targets if graph.has_node(t)]
    #     print(f"Targets in graph: {len(targetes_in_graph)}/{len(targets)}")

    #     # Check reachability
    #     reachability, reachable_count, total_targets = check_target_reachability(graph, sources, targets)
    #     print(f"Target reachability: {reachable_count}/{total_targets} ({reachability*100:.1f}%)")
        
    #     if reachability >= CONNECTIVITY_THRESHOLD:
    #         print(f"Success on attempt {attempt}")

    #         output_dir = processed_directory / pathway / "interactomes"
    #         output_dir.mkdir(parents=True, exist_ok=True)
            
    #         merged_df = pandas.DataFrame()
    #         output_path = output_dir / f"interactome_{args.threshold}.txt"
    #         merged_df.to_csv(output_path, sep="\t", index=False, header=False)
            
    #         break

    #     # Of all original source→target pairs with paths, how many still have paths?
    #     connected_pairs = len(find_connected_sources_targets(sources, targets, graph))
    #     original_pairs = len(find_connected_sources_targets(sources, targets, original_graph))
    #     connectivity = connected_pairs / original_pairs if original_pairs > 0 else 0
    #     print(f"connectivity: {connectivity}")
                
    #     print(f"Below threshold ({CONNECTIVITY_THRESHOLD*100:.0f}%), retrying (attempt {attempt})...")
    #     sampled_interactome = sample_interactome(interactome_df, weight_mapping, args.threshold)
    #     merged_df = pathway_df.merge(
    #         sampled_interactome[["Node1", "Node2"]], 
    #         how="inner", 
    #         on=["Node1", "Node2"]
    #     )
    #     attempt += 1
    
    # Retry loop
    CONNECTIVITY_THRESHOLD = 0.8  # 80%
    attempt = 1
    merged_df = None
    while True:
        sampled_interactome = sample_interactome(interactome_df, weight_mapping, args.threshold)
        
        # print(f"\n--- Attempt {attempt} ---")
        # print(f"Sampled interactome shape: {sampled_interactome.shape}")
        # print(f"Sampled interactome sample rows:\n{sampled_interactome.head(3)}")
        # print(f"Unique Node1 in sample: {sampled_interactome['Node1'].nunique()}")
        # # Check if samples are actually different
        # sample_hash = pandas.util.hash_pandas_object(sampled_interactome, index=True).sum()
        # print(f"Attempt {attempt} - Sample hash: {sample_hash}")
        # print(f"Sample rows 100-105:\n{sampled_interactome.iloc[100:105]}")
        
        median_weight = sampled_interactome["Weight"].median()
        pathway_df["Weight"] = median_weight
        
        original_graph = networkx.from_pandas_edgelist(
            pathway_df, 
            source="Node1", 
            target="Node2", 
            create_using=networkx.DiGraph
        )

        merged_df = attempt_sample(pathway_df, sampled_interactome, sources, targets, original_graph, CONNECTIVITY_THRESHOLD)
        if merged_df is not None:
            break
        attempt += 1

    # Save output
    output_dir = processed_directory / pathway / "interactomes"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"interactome_{args.threshold}.txt"
    merged_df.to_csv(output_path, sep="\t", index=False, header=False)
    print(f"Success on attempt {attempt - 1}")

if __name__ == "__main__":
    main()
