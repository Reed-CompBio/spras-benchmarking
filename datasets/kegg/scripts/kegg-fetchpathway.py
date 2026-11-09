## code taken from benchmarking repo
## datasets/panther_signaling_pathways/scripts/fetch_pathway.py

import argparse
import json
import urllib.parse

from pathlib import Path

from paxtools.fetch import fetch
from paxtools.sif import toSIF

kegg_directory = Path(__file__).parent.resolve() #modified from panther -> kegg


def parser():
    parser = argparse.ArgumentParser(prog="KEGG pathway fetcher")

    parser.add_argument("pathway_name", type=str, nargs="?", default=None)

    return parser

'''
keeping original code format for easier reference
'''
# def main():
#     args = parser().parse_args()
#     all_pathways_df = json.loads((kegg_directory / "intermediate" / "kegg_pathways.json").read_text()) ##in intermediate folder as remnants of og file
#     associated_id = all_pathways_df[urllib.parse.unquote(args.pathway_name)] #unquotes query

#     pathway_data_dir = kegg_directory / "intermediate" / "pathway-pc-data"
#     pathway_data_dir.mkdir(exist_ok=True, parents=True)

#     fetch(
#         kegg_directory / "raw" / "pc-kegg-biopax.owl",
#         pathway_data_dir / Path(args.pathway_name).with_suffix(".owl"),
#         denylist=kegg_directory / "raw" / "blacklist.txt", #remove blacklist names in file
#         uris=[associated_id],
#         absolute=True,
#     )

#     toSIF( #useful post-just having reconx
#         pathway_data_dir / Path(args.pathway_name).with_suffix(".owl"),
#         pathway_data_dir / Path(args.pathway_name).with_suffix(".sif"),
#         # See the paxtools library for information about how these settings were retrieved.
#         # These are directly from PathwayCommons.
#         denylist=str(kegg_directory / "raw" / "blacklist.txt"),
#         chemDb=["chebi"], #for metabolites
#         seqDb=["hgnc"],
#         exclude=["NEIGHBOR_OF"], 
#         extended=True,
#     )

def process_pathway(pathway_name: str, uri: str):
    pathway_data_dir_owl = kegg_directory / "intermediate" / "owl-kegg-pathway-pc-data"
    pathway_data_dir_owl.mkdir(exist_ok=True, parents=True)

    pathway_data_dir_sif = kegg_directory / "intermediate" / "sif-kegg-pathway-pc-data"
    pathway_data_dir_sif.mkdir(exist_ok=True, parents=True)

    #pathway = urllib.parse.quote(pathway_name, safe="-")
    pathway = urllib.parse.quote(pathway_name.replace(" ", "-").replace(",", "-"), safe="-")
    ## trying to get this to work how I want it to

    fetch(
        kegg_directory / "output" / "pc-kegg-biopax.owl",
        pathway_data_dir_owl / Path(pathway).with_suffix(".owl"),
        denylist=kegg_directory / "raw" / "blacklist.txt",
        uris=[uri],
        absolute=True,
    )

    toSIF(
        pathway_data_dir_owl / Path(pathway).with_suffix(".owl"),
        pathway_data_dir_sif / Path(pathway).with_suffix(".sif"),
        denylist=str(kegg_directory / "raw" / "blacklist.txt"),
        chemDb=["chebi"],
        seqDb=["hgnc"],
        exclude=["NEIGHBOR_OF"],
        extended=True,
    )


def main():
    args = parser().parse_args()
    kegg_pathways = json.loads((kegg_directory / "intermediate" / "kegg_pathways.json").read_text())

    if args.pathway_name:
        # single pathway passed via CLI
        name = urllib.parse.unquote(args.pathway_name)
        if name not in kegg_pathways:
            print(f"Pathway '{name}' not found in kegg_pathways.json.")
            return
        process_pathway(name, kegg_pathways[name])
    else:
        # process all pathways in JSON
        print(f"Processing {len(kegg_pathways)} KEGG pathways...")
        for name, uri in kegg_pathways.items():
            print(f"  - {name}")
            process_pathway(name, uri)

if __name__ == "__main__":
    main()