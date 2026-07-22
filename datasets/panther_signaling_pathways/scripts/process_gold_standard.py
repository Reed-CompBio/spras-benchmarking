import pandas as pd
from pathlib import Path

from datasets.panther_signaling_pathways.scripts.parser import parser

panther_directory = Path(__file__).parent.parent.resolve()
intermediate_directory = panther_directory / "intermediate"
raw_directory = panther_directory / "raw"
processed_directory = panther_directory / "processed"

def main():

    intermediate_directory.mkdir(exist_ok=True)
    processed_directory.mkdir(exist_ok=True)
    pathway = parser().parse_args().pathway
    pathway_folder = intermediate_directory / pathway
    processed_gold_standard_folder = processed_directory / pathway
    processed_gold_standard_folder.mkdir(exist_ok=True)


    # read the pathway, and have the use columns command
    # trim the pathawy to the interactome(s)
    # set as the gold standard
    # save

    interactome_file = processed_directory/ "interactome.tsv"

    interactome = pd.read_csv(interactome_file, sep="\t", header=None, names=["Node1", "Node2", "Weight", "Direction"])
    interactome.drop(columns="Weight", inplace=True)
    gold_standard = pd.read_csv(pathway_folder / "pathway.txt", sep="\t", usecols=["Node1", "Node2", "Direction"])

    # Trim gold standard to edges present in the interactome, need to include the direction but can't do this yet because I haven't added the edges into the interactome yet that is in the pathway
    gold_standard = gold_standard.merge(
        interactome[["Node1", "Node2"]],
        on=["Node1", "Node2"],
        how="inner"
    )


if __name__ == "__main__":
    main()