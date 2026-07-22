import argparse
import json
from pathlib import Path

from paxtools.fetch import fetch
from paxtools.sif import toSIF
import urllib.parse

panther_directory = Path(__file__).parent.parent.resolve()


def parser():
    parser = argparse.ArgumentParser(prog="PANTHER pathway fetcher")

    parser.add_argument("pathway_name", type=str)

    return parser


def main():
    args = parser().parse_args()
    curated_pathways_df = json.loads((panther_directory / "intermediate" / "curated_pathways_id_mapping.json").read_text())
    associated_id = curated_pathways_df[urllib.parse.unquote(args.pathway_name)]

    pathway_data_dir = panther_directory / "intermediate" / "pathway-pc-data"
    pathway_data_dir.mkdir(exist_ok=True, parents=True)

    fetch(
        panther_directory / "raw" / "pc-panther-biopax.owl",
        pathway_data_dir / Path(args.pathway_name).with_suffix(".owl"),
        denylist=panther_directory / "raw" / "blacklist.txt",
        uris=[associated_id],
        absolute=True,
    )

    toSIF(
        pathway_data_dir / Path(args.pathway_name).with_suffix(".owl"),
        pathway_data_dir / Path(args.pathway_name).with_suffix(".sif"),
        # See the paxtools library for information about how these settings were retrieved.
        # These are directly from PathwayCommons.
        denylist=str(panther_directory / "raw" / "blacklist.txt"),
        chemDb=["chebi"],
        seqDb=["hgnc"],
        exclude=["NEIGHBOR_OF"],
        extended=True,
    )

if __name__ == "__main__":
    main()
