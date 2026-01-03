import argparse
import networkx
import os
from pathlib import Path
import pandas as pd

dir_path = Path(os.path.dirname(os.path.realpath(__file__)))

def args():
    parser = argparse.ArgumentParser(
                    prog='WikiPathways parse',
                    description='Parses WikiNetworks directed graphs into ENSP-based directed graphs.')
    parser.add_argument('-p', '--pathway', required=True)
    return parser.parse_args()

def transitive_closure(graph: networkx.DiGraph, vertex):
    """Takes a vertex v and returns the graph without v, where neighbors of v are connected transitively."""
    in_neighbours = list(graph.in_edges(vertex))
    out_neighbours = list(graph.out_edges(vertex))

    if len(in_neighbours) == 0 or len(out_neighbours) == 0: return

    # Take a representative edge without the associated edgeID: we want to make sure that the transitive closure of the following edges is the same
    rep_edge = graph.get_edge_data(in_neighbours[0][0], in_neighbours[0][1])
    rep_edge_data = {k: v for k, v in rep_edge.items() if k != 'edgeID'}

    for (u, _) in in_neighbours:
        for (_, w) in out_neighbours:
            assert {k: v for k,v in graph.get_edge_data(u, vertex).items() if k != 'edgeID'} == rep_edge_data
            assert {k: v for k,v in graph.get_edge_data(vertex, w).items() if k != 'edgeID'} == rep_edge_data
            graph.add_edge(u, w, **rep_edge_data)
    graph.remove_node(vertex)

def map_ensembl(graph: networkx.DiGraph):
    """Maps non-Ensembl genes in a network to Ensembl genes. This mutates the input graph."""
    db_nodes = {k: v for k, v in networkx.get_node_attributes(graph, 'database').items() if v != 'Ensembl'}

    # Check that we only have the databases we can map
    valid_databases = set(['Entrez Gene', 'WikiPathways', 'HMDB'])
    assert set(db_nodes.values()).union(valid_databases) == valid_databases, f"Got {set(db_nodes.values()).union(valid_databases)}"

    # WikiPathways genes are usually extraneous genes: we throw them away
    wikipathway_genes = [k for k, v in db_nodes.items() if v == 'WikiPathways']
    graph.remove_nodes_from(wikipathway_genes)

    # HMDB-marked nodes are intermediate products: we can remove these nodes and link incoming nodes to outgoing nodes and vice versa
    # by substituting the node with its transitive closure.
    hmdb_nodes = [k for k, v in db_nodes.items() if v == 'HMDB']
    for node in hmdb_nodes: transitive_closure(graph, node)

    # We map Entrez genes to their Ensembl (ENSG) counterparts.
    entrez_genes = {
        node: data["databaseID"] for node, data in graph.nodes(data=True) if "database" in data and data["database"] == 'Entrez Gene'
    }
    entres_df = pd.read_csv(dir_path / '..' / 'raw' / 'entrez-ensg.tsv', sep='\t', header=None, names=['ENSG', 'Entrez'], dtype=str)
    for node, gene in entrez_genes.items():
        new_id = entres_df[entres_df["Entrez"] == gene]['ENSG'].iloc[0]
        networkx.set_node_attributes(graph, {node: 'Ensembl'}, 'database')
        networkx.set_node_attributes(graph, {node: new_id}, 'databaseID')

    return graph

def flatten_group(graph: networkx.DiGraph):
    pass # TODO

def main():
    pathway = args().pathway
    graph: networkx.DiGraph = networkx.read_graphml(dir_path / ".." / "processed" / f"{pathway}.graphml")
    
    # First, we need to map all non-Ensembl genes to Ensembl genes
    map_ensembl(graph)

    # Then, we take the strongly connected group components and turn them into a proper group
    # representation (such that groups are flattened and not connected to each other, but rather,
    # they form a node layer mimicking that of a hyperedge)
    flatten_group(graph)

    networkx.write_graphml(graph, dir_path / '..' / 'processed' / f"{pathway}-processed.graphml")

if __name__ == "__main__":
    main()
