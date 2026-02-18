"""
Tools for sampling interactomes.
"""

import collections
import networkx
import itertools
import pandas
import random
from typing import OrderedDict, Optional
import os

def count_weights(weights: dict[float, int]) -> OrderedDict[float, int]:
    """
    Returns an ordered map (lowest to highest weight) from the
    weight to the number of elements the weight has.

    The full workflow for this function should be:
    ```python
    count_weights(dict(interactome_df["Weight"].value_counts(sort=False).values))
    ```
    """
    return collections.OrderedDict(sorted({k: int(v) for k, v in weights.items()}.items()))

def find_connected_sources_targets(
        sources: list[str],
        targets: list[str],
        graph: networkx.Graph
) -> list[tuple[str, str]]:
    connections: list[tuple[str, str]] = []
    for source, target in itertools.product(sources, targets):
        if graph.has_node(source) and graph.has_node(target) and networkx.has_path(graph, source, target):
            connections.append((source, target))
    return connections

def attempt_sample(
    pathway_name: str,
    pathway_df: pandas.DataFrame,
    percentage: float,
    weight_mapping: OrderedDict[int, int],
    interactome_df: pandas.DataFrame,
    sources: list[str],
    targets: list[str],
    output_interactome: str | os.PathLike,
    output_gold_standard: str | os.PathLike
) -> Optional[list[tuple[str, str]]]:
    # TODO: generalize to node prizes/actives
    """
    Samples a {pathway_df} (logged as {pathway_name}) along with its backing {interactome_df}
    with a certain {percentage} backed by a {weight_mapping} while preserving some {sources} and {targets},
    outputting to {output_interactome} and {output_gold_standard},
    returning the connections between {sources} and {targets},
    or None if the target percentage failed.
    """
    interactome_df = sample_interactome(interactome_df, weight_mapping, percentage)

    print(f"Merging {pathway_name} with interactome...")
    # While we are merging this graph, we are preparing to compare the connectedness of the prev[ious] and curr[ent] (merged) graph.
    prev_graph = networkx.from_pandas_edgelist(pathway_df, source="Interactor1", target="Interactor2")
    prev_connections = find_connected_sources_targets(sources, targets, prev_graph)

    print("Checking for pathway connectedness...")
    pathway_df = pathway_df.merge(interactome_df, how="inner", on=["Interactor1", "Interactor2"])
    curr_graph = networkx.from_pandas_edgelist(pathway_df, source="Interactor1", target="Interactor2")
    curr_connections = find_connected_sources_targets(sources, targets, curr_graph)

    # We ask that at least `percentage` of the sources and targets are connected with one another.
    connection_percentage = float(len(curr_connections)) / float(len(prev_connections))

    if percentage < connection_percentage:
        print(f"Got {connection_percentage * 100:.1f}% connections above the {percentage * 100:.1f}% threshold.")
        pathway_df.to_csv(output_gold_standard, sep="\t", index=False, header=False)
        interactome_df.to_csv(output_interactome, sep='\t', index=False, header=False)
        return curr_connections
    print(f"Failed {connection_percentage * 100:.1f}% connections below the {percentage * 100:.1f}% threshold.")
    return None

def sample_interactome(
        interactome_df: pandas.DataFrame,
        weight_mapping: OrderedDict[int, int],
        percentage: float
):
    """
    Samples X% of an interactome using its weight_counts dictionary. (See `count_weights` for generating `weight_counts`.)
    """
    if percentage > 1:
        raise RuntimeError(f"Got a percentage above 1 ({percentage})?")
    if percentage == 1:
        return interactome_df
    # Using a list then creating the set is faster because of the sets rather than the gets.
    print("Creating item samples...")
    full_list: list[int] = []
    curr_v = 0
    for k, v in weight_mapping.items():
        full_list.extend(map(lambda x: x + curr_v, random.sample(range(1, v), round(percentage * v))))
        curr_v += v
    full_set = set(full_list)

    print("Sampling interactome...")
    return interactome_df.iloc[list(full_set)]
