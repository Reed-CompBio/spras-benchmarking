from pathlib import Path
from matplotlib import pyplot as plt
import pandas as pd
import networkx as nx
import pickle
import sys


def read_file(path, delimiter):
    try:
        # Read the TSV file
        df = pd.read_csv(path, sep=delimiter)
        return df
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None


def export_graph_to_pickle(graph, filename):
    with open(filename, "wb") as f:
        pickle.dump(graph, f)


def import_graph_from_pickle(filename):
    with open(filename, "rb") as f:
        return pickle.load(f)


def create_interactome(data, columns, threshold=1):
    G = nx.DiGraph()
    for index, row in data.iterrows():
        node1 = row[columns[0]]
        node2 = row[columns[1]]
        experiment = None
        if len(columns) == 3:
            experiment = row[columns[2]]
        else:
            experiment = 1
        if experiment >= threshold:
            if not G.has_node(node1):
                G.add_node(node1)
            if not G.has_node(node2):
                G.add_node(node2)

            if not G.has_edge(node1, node2) or not G.has_edge(node2, node1):
                G.add_edge(node1, node2, weight=experiment)
                G.add_edge(node2, node1, weight=experiment)
    return G


def bfs_k_hop(G: nx.DiGraph, node_set, k):
    k_hop_neighbors_set = set()

    # Iterate over each node in the node set
    for node in node_set:
        # Check if the node exists in the graph
        if node in G:
            # Perform BFS with depth limit k
            bfs_edges = nx.bfs_edges(G, node, depth_limit=k)

            # Extract nodes from the BFS traversal
            bfs_nodes = set([node for edge in bfs_edges for node in edge])
            bfs_nodes.add(node)  # Include the starting node itself

            # Add the BFS nodes to the global set
            k_hop_neighbors_set.update(bfs_nodes)

    return k_hop_neighbors_set


def main():
    pathway_dir = "./pathways"
    stringdb_path = "9606.protein.links.detailed.v12.0.txt"
    # pathways = ["apoptosis", "hedgehog", "jak", "notch", "wnt", "beta3"]
    pathways = ["beta3"]

    interactome_pickle_path = Path("interactome-experimental.pickle")

    # read in files
    df_stringdb = read_file(stringdb_path, " ")

    columns = ["protein1", "protein2", "experimental"]
    threshold = 300
    G: nx.DiGraph = import_graph_from_pickle(interactome_pickle_path)
    # G: nx.DiGraph = create_interactome(df_stringdb, columns, threshold)
    # export_graph_to_pickle(G, f"interactome-experimental-{threshold}.pickle")
    print("Interactome")
    print("nodes:" + str(len(G.nodes())))
    print("edges:" + str(len(G.edges())))
    print()

    for path in pathways:
        file = pathway_dir + "/" + path + "/final.csv"
        df_pathway = read_file(file, ",")

        columns = [
            "string_id1",
            "string_id2",
        ]
        # D: nx.DiGraph = import_graph_from_pickle("pathway.pickle")
        D: nx.DiGraph = create_interactome(df_pathway, columns)
        # export_graph_to_pickle(D, "pathway.pickle")
        print("Pathway: ", path)
        print("nodes:" + str(len(D.nodes())))
        print("edges:" + str(len(D.edges())))
        print()

        print("Pathway and interactome compatibility")
        i = 0
        print(len(D.edges()))
        for edge in D.edges():
            if G.has_edge(edge[0], edge[1]) or G.has_edge(edge[1], edge[0]):
                i += 1

        print("edges from pathway found in interactome: ", i)
        missing = len(D.edges()) - i
        print("missing pathway edges in interactome : ", missing)
        print()

        weakly_connected_components = list(nx.weakly_connected_components(D))

        num_weakly_connected_components = len(weakly_connected_components)

        print("Number of weakly connected components:", num_weakly_connected_components)

        strongly_connected_components = list(nx.strongly_connected_components(D))

        num_strongly_connected_components = len(strongly_connected_components)

        print(
            "Number of strongly connected components:",
            num_strongly_connected_components,
        )
        print()

        # Check if we have disconnected components
        if len(strongly_connected_components) != 1:
            # Create a list of (component, number of edges) tuples
            edge_counts = [
                (component, len(D.subgraph(component).edges()))
                for component in strongly_connected_components
            ]

            # Find the component with the maximum number of edges
            largest_component, max_edges = max(edge_counts, key=lambda x: x[1])
            D = D.subgraph(largest_component)

            print(
                f"Selected component with {max_edges} edges and {len(D.nodes())} nodes."
            )

        plt.figure(figsize=(10, 8))
        pos = nx.spring_layout(D)

        nx.draw(
            D,
            pos,
            with_labels=True,
            node_color="lightblue",
            node_size=700,
            font_size=10,
            font_color="black",
            edge_color="gray",
        )

        plt.title(path)
        output_dir = Path("./images/" + path + ".png")
        plt.savefig(output_dir, format="png", dpi=300)
        plt.show()
        plt.close()


if __name__ == "__main__":
    main()
