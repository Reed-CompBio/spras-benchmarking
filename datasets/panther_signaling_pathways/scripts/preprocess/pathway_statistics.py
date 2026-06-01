"""
Reports on pathway statistics located under `processed`.
"""

from pathlib import Path
import pandas
import urllib.parse
import networkx
import itertools

# from tools.sample.sampling import find_connected_sources_targets

def find_connected_sources_targets(sources: list[str], targets: list[str], graph: networkx.DiGraph) -> list[tuple[str, str]]:
    connections: list[tuple[str, str]] = []
    for source, target in itertools.product(sources, targets):
        if graph.has_node(source) and graph.has_node(target) and networkx.has_path(graph, source, target):
            connections.append((source, target))
    return connections

def convert_undirected_to_directed(df: pandas.DataFrame) -> pandas.DataFrame:
    mask = df["Direction"] == "U"
    new_df = df[mask].copy(deep=True)
    new_df["Interactor1"], new_df["Interactor2"] = new_df["Interactor2"], new_df["Interactor1"]
    new_df["Direction"] = "D"
    df.loc[mask, "Direction"] = "D"
    df = pandas.concat([df, new_df], ignore_index=True)
    return df


current_directory = Path(__file__).parent.resolve()
synthetic_directory = current_directory / ".." / ".."

print(synthetic_directory)

def main():
    data_entries = []

    # We identify pathways by their gold standard edges, since we have a few other files mixed in with `processed`.
    for pathway_folder in (synthetic_directory / "processed" / "pathways").rglob("*/"):
        print(pathway_folder)
        gs_edges_graph = networkx.from_pandas_edgelist(
            convert_undirected_to_directed(
                pandas.read_csv(pathway_folder / "gs_edges.txt", sep="\t", names=["Interactor1", "Interactor2", "Rank", "Direction"])
            ),
            "Interactor1",
            "Interactor2",
            create_using=networkx.DiGraph,
        )
        node_prizes = pandas.read_csv(pathway_folder / "node_prizes.txt", sep="\t")

        sources = list(node_prizes[node_prizes["sources"] == True]["NODEID"])
        targets = list(node_prizes[node_prizes["targets"] == True]["NODEID"])

        n_sources = len(sources)
        n_targets = len(targets)

        connected_sources_targets = find_connected_sources_targets(
            sources,
            targets,
            gs_edges_graph,
        )

        sources_to_targets = (float(n_sources) / float(n_targets)) if n_targets != 0 else float("nan")
        targets_to_sources = (float(n_targets) / float(n_sources)) if n_sources != 0 else float("nan")

        data_entries.append(
            (
                urllib.parse.unquote(pathway_folder.stem),
                n_sources,
                n_targets,
                sources_to_targets,
                targets_to_sources,
                (float(len(connected_sources_targets)) / float(n_sources * n_targets)) if n_sources * n_targets != 0 else 0.0,
            )
        )

    data_df = pandas.DataFrame(
        data_entries,
        columns=("Name", "Sources", "Targets", "Sources/Targets", "Targets/Sources", "Connected Percentage"),
    )
    data_df.to_csv(current_directory / "full_stats.tsv", sep="\t", index=False)

    filtered_df = data_df.loc[data_df["Sources"] != 0].loc[data_df["Targets"] != 0].loc[data_df["Connected Percentage"] != 0]
    filtered_df.reset_index(inplace=True, drop=True)
    print(filtered_df)

    filtered_df.to_csv(current_directory / "filtered_stats.tsv", sep="\t", index=False)



if __name__ == "__main__":
    main()