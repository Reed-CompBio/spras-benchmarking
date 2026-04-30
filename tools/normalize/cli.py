"""
Allows for easy SPRAS post-data processing for any dataset via a CLI callable in Snakemake.
Specifically, this:
- Trims other input files that aren't the interactome.
- Removes duplicate edges from the interactome.

See EGFR for a Snakemake dataset-based example usage of this.
"""

import argparse
from pathlib import Path
import pandas

from tools.normalize.interactome import normalize_interactome
from tools.normalize.trim_input_nodes import trim_input_nodes_file
from tools.normalize.trim_list import trim_node_list_file


def argparser():
    # The language and setup of `deciding_parser` and `parser` is from
    # https://stackoverflow.com/a/70716254/7589775. This allows us to make
    # the other outputs conditionally required.
    deciding_parser = argparse.ArgumentParser(add_help=False)
    deciding_parser.add_argument("--interactome", help="The input interactome file", required=True)
    deciding_parser.add_argument("--interactome-output", help="Where to output the new interactome file.", required=True)
    deciding_parser.add_argument("--input-nodes", help="The input nodes file to trim against.")
    deciding_parser.add_argument("--gold-standard", help="The gold standard file to trim against.")

    deciding_args, _ = deciding_parser.parse_known_args()

    parser = argparse.ArgumentParser(parents=[deciding_parser])
    parser.add_argument("--input-nodes-output", help="Where to output the new input nodes file.", required=deciding_args.input_nodes)
    parser.add_argument("--gold-standard-output", help="Where to output the new gold standard file.", required=deciding_args.gold_standard)

    return parser


def main():
    args = argparser().parse_args()

    print("Normalizing interactome...")
    interactome_df = pandas.read_csv(args.interactome, sep="\t", header=None, names=["Interactor1", "Interactor2", "Weight", "Direction"])
    normalized_interactome, _ = normalize_interactome(interactome_df)
    Path(args.interactome_output).parent.mkdir(parents=True, exist_ok=True)
    normalized_interactome.to_csv(args.interactome_output, sep="\t", index=False)

    if args.input_nodes:
        print("Trimming input nodes...")
        new_input_nodes = trim_input_nodes_file(interactome_df, pandas.read_csv(args.input_nodes, sep="\t"))
        Path(args.input_nodes_output).parent.mkdir(parents=True, exist_ok=True)
        new_input_nodes.to_csv(args.input_nodes_output, sep="\t", index=False)

    if args.gold_standard:
        print("Trimming gold standard nodes list...")
        Path(args.gold_standard_output).parent.mkdir(parents=True, exist_ok=True)
        trim_node_list_file(interactome_df, Path(args.gold_standard), Path(args.gold_standard_output))


if __name__ == "__main__":
    main()
