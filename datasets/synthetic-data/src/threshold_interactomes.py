import pandas as pd
from pathlib import Path
import os

current_directory = Path(os.path.dirname(os.path.realpath(__file__)))
data_directory = current_directory / '..'

if not (data_directory / "processed" / "interactomes" / "uniprot-threshold-interactomes").exists():
    (data_directory / "processed" / "interactomes" / "uniprot-threshold-interactomes").mkdir(exist_ok=True, parents=True)

# get the string -> uniprot accession ID pairings
UniProt_AC = pd.read_csv(data_directory / "raw" / "human-interactome" / "String_to_Uniprot_ids_2025_04_06.tsv", sep='\t', header=0)
one_to_many_dict = UniProt_AC.groupby("From")["Entry"].apply(list).to_dict()

# read in interactome
human_interactome = pd.read_csv(data_directory / "intermediate" / '9606.protein.links.full.v12.0.txt', sep=' ', header = 0)

def get_aliases(protein_id):
    return one_to_many_dict.get(protein_id, [])

human_interactome['protein1_uniprot'] = human_interactome['protein1'].apply(get_aliases)
human_interactome['protein2_uniprot'] = human_interactome['protein2'].apply(get_aliases)

human_interactome = human_interactome.explode('protein1_uniprot').explode('protein2_uniprot')

missing_alias_edges = human_interactome[(human_interactome['protein1_uniprot'].isna()) | (human_interactome['protein2_uniprot'].isna())]

proteins_without_aliases = pd.concat([
    missing_alias_edges.loc[missing_alias_edges['protein1_uniprot'].isna(), 'protein1'],
    missing_alias_edges.loc[missing_alias_edges['protein2_uniprot'].isna(), 'protein2']
], ignore_index=True).drop_duplicates().reset_index(drop=True)
proteins_without_aliases = proteins_without_aliases.to_frame(name="protein")

removed_edges = missing_alias_edges[['protein1', 'protein2']]
removed_edges = removed_edges.drop_duplicates().reset_index(drop=True)

proteins_without_aliases.to_csv(
    data_directory / "processed" / "interactomes" / "uniprot-threshold-interactomes/proteins_missing_aliases.csv",
    sep='\t', index=False, header=True)
removed_edges.to_csv(data_directory / "processed" / "interactomes" / "uniprot-threshold-interactomes/removed_edges.txt", sep='\t', index=False, header=True)

human_interactome = human_interactome.dropna(subset=['protein1_uniprot', 'protein2_uniprot']).reset_index(drop=True)

# threshold the interactomes
thresholds = [1, 100, 200, 300, 400, 500, 600, 700, 800, 900]
for thresh in thresholds:
    thresh_df = human_interactome[human_interactome['experiments'] >= thresh]
    thresh_df = thresh_df[['protein1_uniprot', 'protein2_uniprot', "experiments"]]
    thresh_df["Direction"] = "U"
    thresh_df.columns = ["Node1", "Node2", "Rank", "Direction"]

    # sort by rank, then by (node1 and node2) to ensure deterministic sorting
    thresh_df = thresh_df.sort_values(by=["Rank", "Node1", "Node2"], ascending=True, ignore_index=True)

    # for undirected edges, sort node pairs so that Node1 is always the lesser of the two
    undirected_mask = thresh_df["Direction"] == "U"

    # computes the minimum and maximum node (sorted order) for each row under the mask
    min_nodes = thresh_df.loc[undirected_mask, ["Node1", "Node2"]].min(axis=1)
    max_nodes = thresh_df.loc[undirected_mask, ["Node1", "Node2"]].max(axis=1)

    # assigns the sorted Node1 and Node2 back into the df
    thresh_df.loc[undirected_mask, "Node1"] = min_nodes
    thresh_df.loc[undirected_mask, "Node2"] = max_nodes

    # keep highest rank version of the edge, drop all the others
    thresh_df = thresh_df.sort_values(by=["Node1", "Node2", "Rank"], ascending=[True, True, False], ignore_index=True)
    thresh_df = thresh_df.drop_duplicates(subset=["Node1", "Node2"], keep="first")

    thresh_df.to_csv(
        data_directory / "processed" / "interactomes" / "uniprot-threshold-interactomes" / f"uniprot_human_interactome_{thresh}.txt",
        sep = "\t", header=False, index=False)
