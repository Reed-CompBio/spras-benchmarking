import io
import pandas as pd
from pathlib import Path

from datasets.panther_signaling_pathways.scripts.parser import parser
from tools.trim import trim_data_file

panther_directory = Path(__file__).parent.parent.resolve()
pathway_pc_data_directory = panther_directory / "intermediate" / "pathway-pc-data"
intermediate_directory = panther_directory / "intermediate"

# edge inferred binary relations from Pathway Commons https://www.pathwaycommons.org/pc2/formats
directed = [
    "controls-state-change-of",
    "controls-transport-of",
    "controls-phosphorylation-of",
    "controls-expression-of",
    "catalysis-precedes",
    "consumption-controlled-by",
    "controls-production-of",
    "controls-transport-of-chemical",
    "chemical-affects",
    "used-to-produce",
]

undirected = ["in-complex-with", "interacts-with", "neighbor-of", "reacts-with"]

def process_pathway(pathway_file: Path, folder: Path):
    pathway_file_content = pathway_file.read_text()
    print(pathway_file_content)
    # This file has two csv files stacked on top of each other.
    # This is the header that we are looking for
    needle = "PARTICIPANT\tPARTICIPANT_TYPE\tPARTICIPANT_NAME\tUNIFICATION_XREF\tRELATIONSHIP_XREF"

    edges, nodes = pathway_file_content.split(needle)
    # Re-add the header
    nodes = needle + nodes
    # https://stackoverflow.com/a/65018984/7589775 read the text as a file
    pathway_df = pd.read_csv(io.StringIO(edges), header=0, sep="\t")
    nodes_df = pd.read_csv(io.StringIO(nodes), header=0, sep="\t")

    # First, get the relevant info from the edges
    pathway_df = pathway_df[["PARTICIPANT_A", "INTERACTION_TYPE", "PARTICIPANT_B"]]
    pathway_df.columns = ["NODE1", "INTERACTION_TYPE", "NODE2"]
    # removing ChEBI identifiers: these aren't proteins and we therefore are not interested in them.
    pathway_df = pathway_df[~pathway_df["NODE1"].str.startswith("chebi:")]
    pathway_df = pathway_df[~pathway_df["NODE2"].str.startswith("chebi:")]
    pathway_df["Direction"] = pathway_df["INTERACTION_TYPE"].apply(
        lambda x: "D" if x in directed else ("U" if x in undirected else raise_unknown_direction(x))
    )

    # Do the same for the nodes
    nodes_df = nodes_df[["PARTICIPANT", "PARTICIPANT_NAME", "UNIFICATION_XREF"]]
    nodes_df.columns = ["Gene_name", "UniProt_Entry_Name", "Uniprot_Accession_Number"]
    # removing the chebi: prefix
    nodes_df = nodes_df[~nodes_df["Gene_name"].str.startswith("chebi:")]
    # and remove the uniprot: prefix
    nodes_df["Uniprot_Accession_Number"] = nodes_df["Uniprot_Accession_Number"].str.removeprefix("uniprot:")

    # # Save edges and nodes
    pathway_df.to_csv(folder / "pathway.txt", header=True, index=False, sep="\t")
    nodes_df.to_csv(folder / "pathway_nodes.txt", header=True, index=False, sep="\t")

if __name__ == "__main__":
    pathway = parser().parse_args().pathway
    pathway_file = pathway_pc_data_directory / Path(pathway).with_suffix(".sif")
    pathway_folder = panther_directory / "intermediate" / pathway
    print(pathway_folder)
    pathway_folder.mkdir(parents=True, exist_ok=True)
    process_pathway(pathway_file, pathway_folder)