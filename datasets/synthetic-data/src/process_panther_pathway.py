import io
import pandas as pd
from pathlib import Path
import os

current_directory = Path(os.path.dirname(os.path.realpath(__file__)))

pathways = ["Apoptosis_signaling", "B_cell_activation",
            "Beta3_adrenergic_rec", "Cadherin_signaling",
            "Hedgehog_signaling", "Insulin_IGF",
            "Interleukin_signaling", "Notch_signaling",
            "PDGF_signaling", "Ras", "T_cell_activation",
            "Toll_signaling", "Wnt_signaling", "p38_MAPK",
            "Nicotinic_acetylchol", "Fas_signaling",
            "FGF_signaling", "Interferon_gamma_signaling",
            "JAK_STAT_signaling", "VEGF_signaling"]

data_directory = current_directory / '..' / 'raw' / 'pathway-data'
interactome_folder = current_directory / '..' / 'raw' / 'human-interactome'

def process_pathway(file: Path, folder: Path):
  file_content = file.read_text()
  # This file has two csv files stacked on top of each other. This is the header that we are looking for
  needle = "PARTICIPANT\tPARTICIPANT_TYPE\tPARTICIPANT_NAME\tUNIFICATION_XREF\tRELATIONSHIP_XREF"
  
  edges, nodes = file_content.split(needle)
  # Re-add the header
  nodes = needle + nodes
  # https://stackoverflow.com/a/65018984/7589775 read the text
  # as a file.
  edges_df = pd.read_csv(io.StringIO(edges), header=0, sep='\t')
  nodes_df = pd.read_csv(io.StringIO(nodes), header=0, sep='\t')

  # First, get the relevant info from the edges
  edges_df = edges_df[["PARTICIPANT_A", "INTERACTION_TYPE", "PARTICIPANT_B"]]
  edges_df.columns = ["NODE1", "INTERACTION_TYPE", "NODE2"]
  # removing the chebi: prefix
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
  edges_df.to_csv(folder / 'EDGES.txt', header=True, index=False, sep='\t')
  nodes_df.to_csv(folder / 'NODES.txt', header=True, index=False, sep='\t')

  # Then, we need to get the sources and targets, save them,
  # and mark them with 1.0 prizes:

  # First, for our sources, or transcription factors
  human_tfs = pd.read_csv(interactome_folder / 'Homo_sapiens_TF_Uniprot.txt', sep='\t')
  # TODO
  human_tfs.join(how='inner', other=nodes_df, on='')

  surfaceome = pd.read_csv(interactome_folder / 'Homo_sapiens_surfaceome.txt', sep='\t')

if __name__ == '__main__':
  for pathway in pathways:
    pathway_file = data_directory / pathway / Path(pathway).with_suffix('.txt')
    process_pathway(pathway_file, pathway_file.parent)