import argparse
from pathlib import Path

import pandas
from paxtools.fetch import fetch
from paxtools.sif import toSIF

synthetic_directory = Path(__file__).parent.parent.resolve()


def parser():
    parser = argparse.ArgumentParser(prog="PANTHER pathway fetcher")

    parser.add_argument("pathway_name", type=str)

    return parser


def main():
    args = parser().parse_args()
    curated_pathways_df = pandas.read_csv(synthetic_directory / "intermediate" / "curated_pathways.tsv", sep="\t")
    associated_id = curated_pathways_df.loc[curated_pathways_df["Name"] == args.pathway_name].reset_index(drop=True).loc[0]["ID"]

    pathway_data_dir = synthetic_directory / "intermediate" / "pathway-data"
    pathway_data_dir.mkdir(exist_ok=True, parents=True)

    fetch(
        synthetic_directory / "raw" / "pc-panther-biopax.owl",
        pathway_data_dir / Path(args.pathway_name).with_suffix(".owl"),
        denylist=synthetic_directory / "raw" / "denylist.txt",
        uris=[associated_id],
        absolute=True,
    )

    toSIF(
        pathway_data_dir / Path(args.pathway_name).with_suffix(".owl"),
        pathway_data_dir / Path(args.pathway_name).with_suffix(".sif"),
        # See the paxtools library for information about how these settings were retrieved.
        # These are directly from PathwayCommons.
        denylist=str(synthetic_directory / "raw" / "denylist.txt"),
        chemDb=["chebi"],
        seqDb=["hgnc"],
        exclude=["NEIGHBOR_OF"],
        extended=True,
    )


if __name__ == "__main__":
    main()
