import argparse
import io
import pandas as pd
from pathlib import Path

current_directory = Path(__file__).parent.resolve()

data_directory = current_directory / ".." / "raw" / "pathway-data"
interactome_folder = current_directory / ".." / "raw" / "human-interactome"


def process_pathway(file: Path, folder: Path):
    file_content = file.read_text()
    # This file has two csv files stacked on top of each other.
    # This is the header that we are looking for
    needle = "PARTICIPANT\tPARTICIPANT_TYPE\tPARTICIPANT_NAME\tUNIFICATION_XREF\tRELATIONSHIP_XREF"

    edges, nodes = file_content.split(needle)
    # Re-add the header
    nodes = needle + nodes
    # https://stackoverflow.com/a/65018984/7589775 read the text
    # as a file.
    edges_df = pd.read_csv(io.StringIO(edges), header=0, sep="\t")
    nodes_df = pd.read_csv(io.StringIO(nodes), header=0, sep="\t")

    # First, get the relevant info from the edges
    edges_df = edges_df[["PARTICIPANT_A", "INTERACTION_TYPE", "PARTICIPANT_B"]]
    edges_df.columns = ["NODE1", "INTERACTION_TYPE", "NODE2"]
    # removing ChEBI identifiers: these aren't proteins and we therefore are not interested in them.
    edges_df = edges_df[~edges_df["NODE1"].str.startswith("chebi:")]
    edges_df = edges_df[~edges_df["NODE2"].str.startswith("chebi:")]

    # Do the same for the nodes
    nodes_df = nodes_df[["PARTICIPANT", "UNIFICATION_XREF"]]
    nodes_df.columns = ["NODE", "uniprot"]
    # removing the chebi: prefix
    nodes_df = nodes_df[~nodes_df["NODE"].str.startswith("chebi:")]
    # and remove the uniprot: prefix
    nodes_df["uniprot"] = nodes_df["uniprot"].str.removeprefix("uniprot:")

    # Save edges and nodes
    edges_df.to_csv(folder / "EDGES.txt", header=True, index=False, sep="\t")
    nodes_df.to_csv(folder / "NODES.txt", header=True, index=False, sep="\t")

    # Then, we need to get the sources and targets, save them,
    # and mark them with 1.0 prizes:

    # First, for our targets, or transcription factors
    human_tfs = pd.read_csv(interactome_folder / "Homo_sapiens_TF_Uniprot.tsv", sep="\t")
    human_tfs = nodes_df.merge(human_tfs, how="inner", left_on="uniprot", right_on="UniProtKB-AC")
    human_tfs = human_tfs[["NODE", "uniprot"]]
    human_tfs.to_csv(folder / "TARGETS.txt", sep="\t", index=False)

    # Then, for our receptors. NOTE: we skip the first row since it's empty in the XLSX, so this might break if the surfaceome authors fix this.
    human_receptors = pd.read_excel(interactome_folder / "table_S3_surfaceome.xlsx", sheet_name="in silico surfaceome only", skiprows=1)
    human_receptors = human_receptors[["UniProt accession", "Ensembl gene", "Membranome Almen main-class"]]
    human_receptors = human_receptors[human_receptors["Membranome Almen main-class"] == "Receptors"]
    human_receptors = nodes_df.merge(human_receptors, how="inner", left_on="uniprot", right_on="UniProt accession")
    human_receptors = human_receptors[["NODE", "uniprot"]]
    human_receptors.to_csv(folder / "SOURCES.txt", sep="\t", index=False)

    # Finally, scores
    scores = pd.concat([human_tfs, human_receptors]).drop_duplicates()
    scores["prizes"] = 1
    scores["active"] = "true"
    scores.to_csv(folder / "PRIZES.txt", sep="\t", index=False)


def parser():
    parser = argparse.ArgumentParser(prog="PANTHER pathway parser")

    parser.add_argument("pathway", choices=[file.stem for file in data_directory.iterdir()])

    return parser


if __name__ == "__main__":
    pathway = parser().parse_args().pathway
    pathway_file = data_directory / Path(pathway).with_suffix(".txt")
    intermediate_folder = current_directory / ".." / "intermediate" / pathway
    intermediate_folder.mkdir(parents=True, exist_ok=True)
    process_pathway(pathway_file, intermediate_folder)
