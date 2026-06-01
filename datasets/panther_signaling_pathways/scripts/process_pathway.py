import io
import pandas as pd
from pathlib import Path

from datasets.panther_signaling_pathways.scripts.parser import parser
from tools.trim import trim_data_file

panther_directory = Path(__file__).parent.parent.resolve()
pathway_pc_data_directory = panther_directory / "intermediate" / "pathway-pc-data"
intermediate_directory = panther_directory / "intermediate"
raw_directory = panther_directory / "raw"
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

aliases = pd.read_csv(raw_directory / "9606.protein.aliases.txt", sep="\t")

uniprot_ac_aliases = aliases[aliases["source"] == "UniProt_AC"]
uniprot_id_aliases = aliases[aliases["source"] == "UniProt_ID"]
gene_name_aliases = aliases[aliases["source"] == "UniProt_GN_Name"]

ac_to_ensp = uniprot_ac_aliases.groupby("alias")["#string_protein_id"].apply(list).to_dict()
id_to_ensp = uniprot_id_aliases.groupby("alias")["#string_protein_id"].apply(list).to_dict()
gn_to_ensp = gene_name_aliases.groupby("alias")["#string_protein_id"].apply(list).to_dict()

def process_pathway(pathway_file: Path, pathway_folder: Path):

    pathway_file_content = pathway_file.read_text()
    # This file has two csv files stacked on top of each other.
    # This is the header that we are looking for
    needle = "PARTICIPANT\tPARTICIPANT_TYPE\tPARTICIPANT_NAME\tUNIFICATION_XREF\tRELATIONSHIP_XREF"

    edges, nodes = pathway_file_content.split(needle)
    # Re-add the header
    nodes = needle + nodes
    # https://stackoverflow.com/a/65018984/7589775 read the text as a file
    pathway_df = pd.read_csv(io.StringIO(edges), header=0, sep="\t")
    nodes_df = pd.read_csv(io.StringIO(nodes), header=0, sep="\t")


    # Get relevant info for the nodes
    nodes_df = nodes_df[["PARTICIPANT", "PARTICIPANT_NAME", "UNIFICATION_XREF"]]
    nodes_df.columns = ["UniProt_GN_Name", "UniProt_ID", "UniProt_AC"]
    # removing the chebi: prefix
    nodes_df = nodes_df[~nodes_df["UniProt_GN_Name"].str.startswith("chebi:")]
    # and remove the uniprot: prefix
    nodes_df["UniProt_AC"] = nodes_df["UniProt_AC"].str.removeprefix("uniprot:")
    nodes_df.to_csv(pathway_folder / "raw_pathway_nodes.txt", header=True, index=False, sep="\t")

    # Get the ENSP id names
    nodes_df["ENSP"] = nodes_df[["UniProt_GN_Name", "UniProt_ID", "UniProt_AC"]].apply(
        lambda row: ",".join(set(
            gn_to_ensp.get(row["UniProt_GN_Name"], []) +
            id_to_ensp.get(row["UniProt_ID"], []) +
            ac_to_ensp.get(row["UniProt_AC"], [])
        )),
        axis=1
    )
    nodes_df.to_csv(pathway_folder / "pathway_nodes.txt", header=True, index=False, sep="\t")
    # TODO: add a file that has the ids that are empty
    # Save nodes with no ENSP mapping for inspection
    unmatched = nodes_df[nodes_df["ENSP"] == ""]
    if not unmatched.empty:
        unmatched.to_csv(pathway_folder / "unmatched_nodes.txt", header=True, index=False, sep="\t")

    # Get the relevant info from the edges
    pathway_df = pathway_df[["PARTICIPANT_A", "INTERACTION_TYPE", "PARTICIPANT_B"]]
    pathway_df.columns = ["Node1", "INTERACTION_TYPE", "Node2"]
    # removing ChEBI identifiers: these aren't proteins and we therefore are not interested in them.
    pathway_df = pathway_df[~pathway_df["Node1"].str.startswith("chebi:")]
    pathway_df = pathway_df[~pathway_df["Node2"].str.startswith("chebi:")]
    pathway_df.to_csv(pathway_folder / "raw_pathway.txt", header=True, index=False, sep="\t")
    pathway_df["Direction"] = pathway_df["INTERACTION_TYPE"].apply(
        lambda x: "D" if x in directed else ("U" if x in undirected else raise_unknown_direction(x))
    )
    # TODO: convert UniProt_GN_Name to ENSP for the edges
    # - if there is more than one ENSP for an UniProt_GN_Name id, then add edges to have all the versions of the ENSP conversion
    
    # Convert UniProt_GN_Name to ENSP for the edges
    # Map each node to its list of ENSPs (empty list if no mapping)
    pathway_df["Node1"] = pathway_df["Node1"].map(lambda n: gn_to_ensp.get(n, []))
    pathway_df["Node2"] = pathway_df["Node2"].map(lambda n: gn_to_ensp.get(n, []))
    print(pathway_df)

    # Drop edges where either node has no ENSP mapping
    pathway_df = pathway_df[
        pathway_df["Node1"].map(len) > 0
    ]
    pathway_df = pathway_df[
        pathway_df["Node2"].map(len) > 0
    ]

    # Explode both columns to create all combinations
    pathway_df = pathway_df.explode("Node1").explode("Node2").reset_index(drop=True)

    pathway_df.to_csv(pathway_folder / "pathway.txt", header=True, index=False, sep="\t")

    

   

if __name__ == "__main__":
    pathway = parser().parse_args().pathway
    pathway_file = pathway_pc_data_directory / Path(pathway).with_suffix(".sif")
    pathway_folder = panther_directory / "intermediate" / pathway
    pathway_folder.mkdir(parents=True, exist_ok=True)
    process_pathway(pathway_file, pathway_folder)