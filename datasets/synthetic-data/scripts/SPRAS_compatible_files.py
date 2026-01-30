import os
import pandas as pd
from pathlib import Path
import sys

current_directory = Path(os.path.dirname(os.path.realpath(__file__)))

spras_compatible_dir = Path(current_directory, "..", "processed")
directory = Path(current_directory, "..", "intermediate")

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
    "consumption-controled-by"
]

undirected = [
    "in-complex-with",
    "interacts-with",
    "neighbor-of",
    "reacts-with"
]

def raise_unknown_direction(dir: str):
    raise ValueError(f"Unknown direction {dir}")

if __name__ == '__main__':
    if not os.path.exists(spras_compatible_dir):
        os.makedirs(spras_compatible_dir)

    pathway = sys.argv[1]
    pathway_folder = directory / pathway

    # Create the output folder "uniprot" within the pathway directory
    out_folder = os.path.join(spras_compatible_dir, pathway)
    os.makedirs(out_folder, exist_ok=True)

    nodes_file = os.path.join(pathway_folder, "NODES.txt")
    nodes_df = pd.read_csv(nodes_file, sep="\t")

    # a dictionary mapping gene -> Uniprot accession ID
    gene_to_uniprot = pd.Series(
        nodes_df['uniprot'].values,
        index=nodes_df['NODE']).to_dict()

    # nodes
    nodes_uniprot = nodes_df[['uniprot']]
    nodes_uniprot.to_csv(
        os.path.join(out_folder, f"{pathway}_gs_nodes.txt"),
        sep="\t", index=False, header=False)

    # edges
    edges_file = os.path.join(pathway_folder, "EDGES.txt")
    edges_df = pd.read_csv(edges_file, sep="\t", header=0)
    edges_df['NODE1'] = edges_df['NODE1'].map(gene_to_uniprot)
    edges_df['NODE2'] = edges_df['NODE2'].map(gene_to_uniprot)
    edges_df['Rank'] = 1
    edges_df["Direction"] = edges_df["INTERACTION_TYPE"].apply(
        lambda x: "D" if x in directed else
            ("U" if x in undirected else raise_unknown_direction(x))
    )
    edges_df = edges_df.drop(columns = "INTERACTION_TYPE")

    # remove duplicate rows
    # sort by (node1 and node2) to ensure deterministic sorting
    edges_df = edges_df.sort_values(
        by=["NODE1", "NODE2"], ascending=True, ignore_index=True)
    undirected_mask = edges_df["Direction"] == "U"
    min_nodes = edges_df.loc[undirected_mask, ["NODE1", "NODE2"]].min(axis=1)
    max_nodes = edges_df.loc[undirected_mask, ["NODE1", "NODE2"]].max(axis=1)
    edges_df.loc[undirected_mask, "NODE1"] = min_nodes
    edges_df.loc[undirected_mask, "NODE2"] = max_nodes

    # keep 1 directed and 1 undirected edge if both exist
    # since rank is 1, we don't need to sort by rank.
    edges_df = edges_df.sort_values(by=["NODE1", "NODE2", "Direction"],
                                    ascending=True, ignore_index=True)
    edges_df = edges_df.drop_duplicates(keep="first", ignore_index=True)

    edges_df.to_csv(
        os.path.join(out_folder, f"{pathway}_gs_edges.txt"),
        sep="\t", index=False, header=False)

    # prizes, targets, sources
    prizes_file = os.path.join(pathway_folder, "PRIZES.txt")
    prizes_df = pd.read_csv(prizes_file, sep="\t")
    prizes_uniprot = prizes_df[['uniprot', 'prizes', 'active']]

    target_file = os.path.join(pathway_folder, "TARGETS.txt")
    target_df = pd.read_csv(target_file, sep="\t")
    target_uniprot = target_df[['uniprot']]

    source_file = os.path.join(pathway_folder, "SOURCES.txt")
    source_df = pd.read_csv(source_file, sep="\t")
    source_uniprot = source_df[['uniprot']]

    # final resulting df combining all the sources, targets, and prizes
    prizes_df['sources'] = prizes_df['uniprot'].isin(source_df['uniprot'])
    prizes_df['targets'] = prizes_df['uniprot'].isin(target_df['uniprot'])
    prizes_df['dummy'] = ""
    prizes_df.rename(columns={'uniprot': 'NODEID', 'prizes': 'prize'}, inplace=True)
    result_df = prizes_df[['NODEID', 'prize', 'sources', 'targets', 'active', 'dummy']]
    result_df.to_csv(
        os.path.join(out_folder, f"{pathway}_node_prizes.txt"),
        sep="\t", index=False, header=True)
