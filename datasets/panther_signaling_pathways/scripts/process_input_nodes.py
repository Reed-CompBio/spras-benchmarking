# convert to SPRAS format for both files?

import pandas as pd
import pandas as pd
from pathlib import Path
from datasets.panther_signaling_pathways.scripts.parser import parser

panther_directory = Path(__file__).parent.parent.resolve()
intermediate_directory = panther_directory / "intermediate"
raw_directory = panther_directory / "raw"

def process_sources():
    None
    # get the surfacaomes
    # convert all the ids to XXX
    # trim the surfacaomes to what is in the interactome (do this for each of the thresholded interactomes as well)
    # save as an intermediate


def process_targets():
    None
    # get the TFs
    # convert all the ids to XXX
    # trim the tfs to what is in the interactome (do this for each of the thresholded interactomes as well)
    # save as an intermediate


def main():
    pathway = parser().parse_args().pathway

    pathway_folder = intermediate_directory / pathway
    pathway_file = pathway_folder / "pathway.txt"
    pathway_nodes_file = pathway_folder / "pathway_nodes.txt"

    targets = pd.read_csv(raw_directory / "Homo_sapiens_TF.tsv", sep="\t")
    print(targets.head())
    sources = pd.read_excel(raw_directory / "table_S3_surfaceome.xlsx", sheet_name="in silico surfaceome only", skiprows=1)
    print(sources.head())

if __name__ == "__main__":
    main()