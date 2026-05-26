"""
Reports on pathway statistics located under `processed`.
"""

from pathlib import Path
import networkx
import pandas
import urllib.parse

from tools.sample import find_connected_sources_targets


# From SPRAS. TODO: import once SPRAS uses pixi
def convert_undirected_to_directed(df: pandas.DataFrame) -> pandas.DataFrame:
    mask = df["Direction"] == "U"
    new_df = df[mask].copy(deep=True)
    new_df["Interactor1"], new_df["Interactor2"] = new_df["Interactor2"], new_df["Interactor1"]
    new_df["Direction"] = "D"
    df.loc[mask, "Direction"] = "D"
    df = pandas.concat([df, new_df], ignore_index=True)
    return df


current_directory = Path(__file__).parent.resolve()
synthetic_directory = current_directory / ".."


def main():
    data_entries = []

    # We identify pathways by their gold standard edges, since we have a few other files mixed in with `processed`.
    for pathway_folder in (synthetic_directory / "processed" / "pathways").rglob("*/"):
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

        connected_sources_targets = find_connected_sources_targets(
            sources,
            targets,
            gs_edges_graph,
        )
        data_entries.append(
            (
                urllib.parse.unquote(pathway_folder.stem),
                len(sources),
                len(targets),
                (float(len(connected_sources_targets)) / float(len(sources) * len(targets))) if len(sources) * len(targets) != 0 else 0.0,
            )
        )

    data_df = pandas.DataFrame(data_entries, columns=("Name", "Sources", "Targets", "Connected Percentage"))
    data_df.to_csv(current_directory / "full_stats.tsv", sep="\t", index=False)

    filtered_df = data_df.loc[data_df["Sources"] != 0].loc[data_df["Targets"] != 0].loc[data_df["Connected Percentage"] != 0]
    print(filtered_df)


if __name__ == "__main__":
    main()
