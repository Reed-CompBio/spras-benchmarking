import argparse
import pandas
from pathlib import Path
import collections
from typing import OrderedDict, NamedTuple
from tools.sample import attempt_sample

current_directory = Path(__file__).parent.resolve()


# From SPRAS. TODO: import once SPRAS uses pixi
def convert_undirected_to_directed(df: pandas.DataFrame) -> pandas.DataFrame:
    mask = df["Direction"] == "U"
    new_df = df[mask].copy(deep=True)
    new_df["Interactor1"], new_df["Interactor2"] = new_df["Interactor2"], new_df["Interactor1"]
    new_df["Direction"] = "D"
    df.loc[mask, "Direction"] = "D"
    df = pandas.concat([df, new_df], ignore_index=True)
    return df


def parser():
    parser = argparse.ArgumentParser(prog="PANTHER pathway parser")

    parser.add_argument("pathway", choices=[file.stem for file in (current_directory / ".." / "raw" / "pathway-data").iterdir()])

    return parser


def count_weights() -> OrderedDict[int, int]:
    """Returns an ordered map (lowest to highest weight) from the weight to the number of elements the weight has"""
    weight_counts = pandas.read_csv(current_directory / ".." / "processed" / "weight-counts.tsv", sep="\t")
    return collections.OrderedDict(sorted({int(k * 1000): int(v) for k, v in dict(weight_counts.values).items()}.items()))


def read_pathway(pathway_name: str) -> pandas.DataFrame:
    """
    Returns the directed-only pathway from a pathway name,
    with columns Interactor1 -> Interactor2.
    """
    pathway_df = pandas.read_csv(
        current_directory / ".." / "processed" / pathway_name / f"{pathway_name}_gs_edges.txt",
        sep="\t",
        names=["Interactor1", "Interactor2", "Weight", "Direction"],
    )
    # We consider an undirected edge to be two directed edges
    pathway_df = convert_undirected_to_directed(pathway_df)
    return pathway_df[["Interactor1", "Interactor2"]]


class SourcesTargets(NamedTuple):
    sources: list[str]
    targets: list[str]


def sources_and_targets(pathway_name: str) -> SourcesTargets:
    """
    Returns the sources and targets associated with a particular pathway
    """
    nodes_df = pandas.read_csv(
        current_directory / ".." / "processed" / pathway_name / f"{pathway_name}_node_prizes.txt", sep="\t", usecols=["NODEID", "sources", "targets"]
    )
    sources: list[str] = list(nodes_df[nodes_df["sources"] is True]["NODEID"])
    targets: list[str] = list(nodes_df[nodes_df["targets"] is True]["NODEID"])

    return SourcesTargets(sources, targets)


def main():
    pathway_name = parser().parse_args().pathway
    print("Reading interactome...")
    interactome_df = pandas.read_csv(
        current_directory / ".." / "processed" / "interactome.tsv",
        header=None,
        sep="\t",
        names=["Interactor1", "Interactor2", "Weight", "Direction"],
        usecols=[0, 1],
    )

    # For performance reasons (groupby is quite slow), we sample in the interactome using the pre-computed weight-counts.tsv file
    weight_mapping = count_weights()

    # Get information about the pathway
    pathway_df = read_pathway(pathway_name)
    sources, targets = sources_and_targets(pathway_name)

    # TODO: isolate percentage constant (this currently builds up 0%, 10%, ..., 100%)
    for percentage in map(lambda x: (x + 1) / 10, range(10)):
        print(f"Sampling with {percentage * 100:.1f}% of edges...")
        attempt_number = 1
        while attempt_sample(
                pathway_name, pathway_df, percentage,
                weight_mapping, interactome_df, sources, targets,
                # TODO: save attempt number
                # TODO: save sources & targets using trim.py
                output_interactome=current_directory / '..' / 'thresholded' / str(percentage) / pathway_name / 'interactome.txt',
                output_gold_standard=current_directory / '..' / 'thresholded' / str(percentage) / pathway_name / 'gold_standard_edges.txt') is None:
            attempt_number += 1
            print(f"Attempt number {attempt_number}")


if __name__ == "__main__":
    main()
