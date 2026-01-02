import argparse
import networkx
import os
from pathlib import Path
import pandas as pd
import copy
import itertools

dir_path = Path(os.path.dirname(os.path.realpath(__file__)))


def args():
    parser = argparse.ArgumentParser(prog="WikiPathways parse", description="Parses WikiNetworks directed graphs into ENSP-based directed graphs.")
    parser.add_argument("-p", "--pathway", required=True)
    return parser.parse_args()


def transitive_closure(graph: networkx.DiGraph, vertex):
    """Takes a vertex v and returns the graph without v, where neighbors of v are connected transitively."""
    in_neighbours = list(graph.in_edges(vertex))
    out_neighbours = list(graph.out_edges(vertex))

    if len(in_neighbours) == 0 or len(out_neighbours) == 0:
        return

    # Take a representative edge without the associated edgeID: we want to make sure that the transitive closure of the following edges is the same
    rep_edge = graph.get_edge_data(in_neighbours[0][0], in_neighbours[0][1])
    rep_edge_data = {k: v for k, v in rep_edge.items() if k != "edgeID"}

    for u, _ in in_neighbours:
        for _, w in out_neighbours:
            assert {k: v for k, v in graph.get_edge_data(u, vertex).items() if k != "edgeID"} == rep_edge_data
            assert {k: v for k, v in graph.get_edge_data(vertex, w).items() if k != "edgeID"} == rep_edge_data
            graph.add_edge(u, w, **rep_edge_data)
    graph.remove_node(vertex)


def map_ensembl(graph: networkx.DiGraph):
    """Maps non-Ensembl genes in a network to Ensembl genes. This mutates the input graph."""
    db_nodes = {k: v for k, v in networkx.get_node_attributes(graph, "database").items() if v != "Ensembl"}

    # Check that we only have the databases we can map
    valid_databases = set(["Entrez Gene", "WikiPathways", "HMDB"])
    assert set(db_nodes.values()).union(valid_databases) == valid_databases, f"Got {set(db_nodes.values()).union(valid_databases)}"

    # WikiPathways genes are usually extraneous genes: we throw them away
    wikipathway_genes = [k for k, v in db_nodes.items() if v == "WikiPathways"]
    graph.remove_nodes_from(wikipathway_genes)

    # HMDB-marked nodes are secondary messangers: we can remove these nodes and link incoming nodes to outgoing nodes and vice versa
    # by substituting the node with its transitive closure.
    hmdb_nodes = [k for k, v in db_nodes.items() if v == "HMDB"]
    for node in hmdb_nodes:
        transitive_closure(graph, node)

    # We map Entrez genes to their Ensembl (ENSG) counterparts.
    entrez_genes = {node: data["databaseID"] for node, data in graph.nodes(data=True) if "database" in data and data["database"] == "Entrez Gene"}
    # TODO: generalize
    entres_df = pd.read_csv(dir_path / ".." / "raw" / "entrez-ensg.tsv", sep="\t", header=None, names=["ENSG", "Entrez"], dtype=str)
    for node, gene in entrez_genes.items():
        new_id = entres_df[entres_df["Entrez"] == gene]["ENSG"].iloc[0]
        networkx.set_node_attributes(graph, {node: "Ensembl"}, "database")
        networkx.set_node_attributes(graph, {node: new_id}, "databaseID")

    # Drop all nodes with no associated database
    no_db_nodes = set(graph.nodes()) - set(networkx.get_node_attributes(graph, "database").keys())
    for node in no_db_nodes:
        graph.remove_node(node)

    # Normalize node data to just be ENSG: (ensg id)
    for node, data in graph.nodes(data=True):
        databaseID = data["databaseID"]
        data.clear()
        data["ENSG"] = databaseID

    return graph


def map_ensg(graph: networkx.DiGraph):
    ensg_ensp_df = pd.read_csv(dir_path / ".." / ".." / ".." / "databases" / "string" / "9606.protein.aliases.v12.0.txt", sep="\t")
    ensg_ensp_df = ensg_ensp_df[ensg_ensp_df["source"] == "Ensembl_gene"].drop(columns=["source"])

    # We cache the original nodes list since this will be changing over time as we go through the iteration.
    for node, data in list(graph.nodes(data=True)):
        ensg = data["ENSG"]
        data.clear()
        ensps = ensg_ensp_df[ensg_ensp_df["alias"] == ensg]["#string_protein_id"]

        # We have our ENSP IDs: lets delete the original node and make our connected group
        in_edges = ((u, v, data) for u, v, data in graph.in_edges(node, data=True))
        out_edges = ((u, v, data) for u, v, data in graph.out_edges(node, data=True))
        graph.remove_node(node)

        # The ENSP nodes form a strongly connected subgraph
        ensp_nodes = []
        for ensp in ensps:
            ensp = ensp[len("9606.") :]
            graph.add_node(node, ENSP=ensp)
            for other_node in ensp_nodes:
                graph.add_edge(node, other_node, arrow="group")
                graph.add_edge(other_node, node, arrow="group")
            ensp_nodes.append(node)

        # and we re-add the other edges
        for u, _, data in in_edges:
            for ensp_node in ensp_nodes:
                graph.add_edge(u, ensp_node, **data)
        for _, v, data in out_edges:
            for ensp_node in ensp_nodes:
                graph.add_edge(ensp_node, v, **data)


def flatten_group(graph: networkx.DiGraph):
    # We create a group_view, such that only group edges are included
    group_view = networkx.subgraph_view(graph, filter_edge=lambda u, v: graph[u][v]["arrow"] == "group")
    group_view = networkx.subgraph_view(group_view, filter_node=lambda u: networkx.degree(graph, u) != 0)

    # As a quick implementation detail, WikiNetworks always gives bi-directionality to group edges, so we know that
    # attracting components are connected components
    assert networkx.number_connected_components(group_view.to_undirected(as_view=True)) == networkx.number_attracting_components(group_view)
    # There are some grouped nodes that don't form fully strongly connected components: to fix this, we make all attracted components
    # strongly connected components
    for component in networkx.attracting_components(group_view):
        for u in component:
            for v in component:
                if u == v:
                    continue
                if group_view.has_edge(u, v):
                    continue
                graph.add_edge(u, v, arrow="group")

    # Now that we're done with our mutations on the underlying graph, we copy and freeze this as we only care about the static list of nodes.
    # (and mark this as as_view to prevent accidental mutation of group_view)
    group_view = copy.deepcopy(group_view.copy()).copy(as_view=True)

    # Finally, we make the hypergraph-edge-like representations of our components.
    for component in networkx.attracting_components(group_view):
        # Note that we want the in_edges of the _graph_ itself
        in_edges_all = list(itertools.chain(*(list(graph.in_edges(node, data=True)) for node in component)))
        out_edges_all = list(itertools.chain(*(list(graph.out_edges(node, data=True)) for node in component)))
        in_edges = list((u, v, data) for u, v, data in in_edges_all if data["arrow"] != "group")
        out_edges = list((u, v, data) for u, v, data in out_edges_all if data["arrow"] != "group")

        # First, we remove every single in_edges_all and out_edges_all
        for u, v, _data in in_edges_all:
            if graph.has_edge(u, v):
                graph.remove_edge(u, v)
        for u, v, _data in out_edges_all:
            if graph.has_edge(u, v):
                graph.remove_edge(u, v)

        in_nodes = list(((u, data) for u, _, data in in_edges))
        out_nodes = list(((v, data) for _, v, data in out_edges))

        # Then re-associate the in and out edges with our new unconnected group nodes.
        for node in component:
            for u, data in in_nodes:
                graph.add_edge(u, node, **data)
            for v, data in out_nodes:
                graph.add_edge(node, v, **data)


def main():
    pathway = args().pathway
    graph: networkx.DiGraph = networkx.read_graphml(dir_path / ".." / "processed" / f"{pathway}.graphml")

    # We remove self loops generated as artifacts from associated annotations
    graph.remove_edges_from(networkx.selfloop_edges(graph))

    # We need to map all non-Ensembl genes to Ensembl genes
    map_ensembl(graph)
    # and map ENSG to ENSP: we form groups for one-to-many ENSG-ENSP connections.
    map_ensg(graph)

    # We take the connected group components and turn them into a proper group
    # representation (such that groups are flattened and not connected to each other, but rather,
    # they form a node layer mimicking that of a hyperedge)
    flatten_group(graph)

    # (with temporary GraphML output)
    networkx.write_graphml(graph, dir_path / ".." / "processed" / f"{pathway}-processed.graphml")

    # and save our new gold standard file.
    graph = networkx.relabel_nodes(graph, {node: data["ENSP"] for node, data in graph.nodes(data=True)})
    networkx.set_edge_attributes(graph, "D", "Direction")
    networkx.write_edgelist(graph, dir_path / ".." / "processed" / f"gs-{pathway}-unprocessed.tsv", delimiter="\t", data=["Direction"])


if __name__ == "__main__":
    main()
