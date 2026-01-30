import argparse
import itertools
from pathlib import Path
import random
import networkx
import uuid
import pandas

def random_id() -> str:
    return uuid.uuid4().hex

def assign_ids(graph: networkx.DiGraph) -> networkx.DiGraph:
    """Assigns new IDs to a graph based on `random_id`"""
    mapping = {node: random_id() for node in graph}
    return networkx.relabel_nodes(graph, mapping)

def gnp_noise(graph: networkx.DiGraph, p: float):
    """
    The mutative equivalent to networkx.gnp_random_graph,
    whose original implementation does not consume a graph.
    """
    for e in itertools.permutations(graph.nodes, 2):
        if random.random() < p:
            graph.add_edge(*e)

def generate_parser():
    parser = argparse.ArgumentParser(prog='Pathway generator')
    parser.add_argument("--path-count", type=int, default=10, help="The number of paths, whose starts and ends are marked as sources and targets.")
    parser.add_argument("--path-length", type=int, default=7, help="The length of every path from --path-count.")

    parser.add_argument("--sources-output", type=str, default="sources.txt")
    parser.add_argument("--targets-output", type=str, default="targets.txt")

    parser.add_argument("--gold-standard-noise", type=float, default=0.03,
                        help="The probability that edges in the gold standard are connected to each other.")
    parser.add_argument("--gold-standard-output", type=str, default="gold-standard.tsv")

    parser.add_argument("--interactome-extra-nodes", type=int, default=400)
    parser.add_argument("--interactome-noise", type=float, default=0.01,
                        help="The probability that edges in the larger interactome are connected to each other.")
    parser.add_argument("--interactome-output", type=str, default="interactome.tsv")
    return parser

def main():
    args = generate_parser().parse_args()

    graph = networkx.DiGraph()
    sources: list[str] = []
    targets: list[str] = []

    # Add the path graphs to form the base of the pathway, while getting sources and targets as well.
    for _ in range(args.path_count):
        path_graph = networkx.path_graph(args.path_length, create_using=networkx.DiGraph())
        path_graph = assign_ids(path_graph)

        topological_sort = list(networkx.topological_sort(path_graph))
        first_node, last_node = (topological_sort[0], topological_sort[-1])
        sources.append(first_node)
        targets.append(last_node)

        graph = networkx.union(graph, path_graph)

    Path(args.sources_output).write_text("\n".join(sources))
    Path(args.targets_output).write_text("\n".join(targets))

    # Then, we'll add some noise: this will be our gold standard.
    gnp_noise(graph, args.gold_standard_noise)
    gold_standard = pandas.DataFrame(((a, b) for a, b, _data in networkx.to_edgelist(graph)), columns=["Source", "Target"])
    # We make the gold standard output a little annoying to force some post-processing with pandas.
    gold_standard.insert(1, "Interaction-Type", "pp")
    gold_standard.to_csv(args.gold_standard_output, index=False, sep='\t')

    # and we'll follow along similarly to above to build our interactome.
    graph.add_nodes_from((random_id() for _ in range(args.interactome_extra_nodes)))
    gnp_noise(graph, args.interactome_noise)
    interactome = pandas.DataFrame(((a, b) for a, b, _data in networkx.to_edgelist(graph)), columns=["Source", "Target"])
    interactome.to_csv(args.interactome_output, index=False, sep='\t')

if __name__ == "__main__":
    main()
