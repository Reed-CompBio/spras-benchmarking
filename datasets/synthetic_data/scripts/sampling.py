import pandas
from pathlib import Path
import collections
from typing import OrderedDict, NamedTuple

import urllib.parse
from tools.sample import attempt_sample
from tools.trim import trim_data_file
from datasets.synthetic_data.scripts.util.parser import parser
import random

synthetic_directory = Path(__file__).parent.parent.resolve()


# From SPRAS. TODO: import once SPRAS uses pixi
def convert_undirected_to_directed(df: pandas.DataFrame) -> pandas.DataFrame:
    mask = df["Direction"] == "U"
    new_df = df[mask].copy(deep=True)
    new_df["Interactor1"], new_df["Interactor2"] = new_df["Interactor2"], new_df["Interactor1"]
    new_df["Direction"] = "D"
    df.loc[mask, "Direction"] = "D"
    df = pandas.concat([df, new_df], ignore_index=True)
    return df


def count_weights() -> OrderedDict[int, int]:
    """Returns an ordered map (lowest to highest weight) from the weight to the number of elements the weight has"""
    weight_counts = pandas.read_csv(synthetic_directory / "processed" / "weight-counts.tsv", sep="\t")
    return collections.OrderedDict(sorted({int(k * 1000): int(v) for k, v in dict(weight_counts.values).items()}.items()))


def read_pathway(pathway_name: str) -> pandas.DataFrame:
    """
    Returns the directed-coerced pathway from a pathway name,
    with columns Interactor1 -> Interactor2.
    """
    pathway_df = pandas.read_csv(
        synthetic_directory / "processed" / "pathways" / pathway_name / "gs_edges.txt",
        sep="\t",
        names=["Interactor1", "Interactor2", "Weight", "Direction"],
    )
    # We consider an undirected edge to be two directed edges
    pathway_df = convert_undirected_to_directed(pathway_df)
    return pathway_df[["Interactor1", "Interactor2"]]


class SourcesTargets(NamedTuple):
    sources: list[str]
    targets: list[str]


def get_node_data(pathway_name: str) -> pandas.DataFrame:
    return pandas.read_csv(
        synthetic_directory / "processed" / "pathways" / pathway_name / "node_prizes.txt", sep="\t", usecols=["NODEID", "sources", "targets"]
    )


def sources_and_targets(pathway_node_prizes_df: pandas.DataFrame) -> SourcesTargets:
    """
    Returns the sources and targets associated with a particular pathway
    """
    sources: list[str] = list(pathway_node_prizes_df[pathway_node_prizes_df["sources"] == True]["NODEID"])
    targets: list[str] = list(pathway_node_prizes_df[pathway_node_prizes_df["targets"] == True]["NODEID"])

    return SourcesTargets(sources, targets)


def main():
    arg_parser = parser()
    arg_parser.add_argument("--seed", help="The randomness seed to use", type=int, required=False)
    arg_parser.add_argument("--amount", help="The amount of thresholds to use", type=int, default=10)
    arg_parser.add_argument(
        "--percentage_thresholding_multiplier",
        help="The percentage multiplier to threshold by, " + "to unlink the sampling percentage to the actual required percentage of connections",
        type=float,
        default=1.0,
    )

    args = arg_parser.parse_args()
    pathway_location = args.pathway
    pathway_name = urllib.parse.unquote(pathway_location)
    if args.seed is not None:
        random.seed(args.seed)

    print("Reading interactome...")
    interactome_df = pandas.read_csv(
        synthetic_directory / "processed" / "interactome.tsv",
        header=None,
        sep="\t",
        names=["Interactor1", "Interactor2", "Weight", "Direction"],
        usecols=[0, 1],
    )

    # For performance reasons (groupby is quite slow), we sample in the interactome using the pre-computed weight-counts.tsv file
    weight_mapping = count_weights()

    # Get information about the pathway
    pathway_df = read_pathway(pathway_location)
    node_data_df = get_node_data(pathway_location)
    sources, targets = sources_and_targets(node_data_df)

    percentages = list(map(lambda x: (x + 1) / args.amount, range(args.amount)))
    for percentage_to_sample in percentages:
        percentage_to_require = percentage_to_sample * args.percentage_thresholding_multiplier

        output_directory = synthetic_directory / "thresholded" / str(percentage_to_sample) / pathway_location
        output_directory.mkdir(exist_ok=True, parents=True)
        output_interactome = output_directory / "interactome.txt"
        output_gold_standard = output_directory / "gold_standard_edges.txt"

        print(f"Sampling with {percentage_to_sample * 100:.1f}% of edges...")
        attempt_number = 1
        while (
            attempt_sample(
                pathway_name,
                pathway_df,
                percentage_to_sample,
                percentage_to_require,
                weight_mapping,
                interactome_df,
                sources,
                targets,
                output_interactome=output_interactome,
                output_gold_standard=output_gold_standard,
            )
            is None
        ):
            attempt_number += 1
            print(f"Attempt number {attempt_number}")

        # We're done sampling:
        (output_directory / "attempt-number.txt").write_text(str(attempt_number))
        # we need to trim our data file as well. We do this already in process_panther_pathway, though.
        trim_data_file(data_df=node_data_df, gold_standard_df=pathway_df).to_csv(output_directory / "node_prizes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
