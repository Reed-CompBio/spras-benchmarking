import pandas as pd
from pathlib import Path
from datasets.synthetic_data.scripts.util.parser import parser
from tools.trim import trim_data_file

synthetic_directory = Path(__file__).parent.parent.resolve()

processed_directory = synthetic_directory / "processed"
intermediate_directory = synthetic_directory / "intermediate"

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
    "consumption-controled-by",
]

undirected = ["in-complex-with", "interacts-with", "neighbor-of", "reacts-with"]


def raise_unknown_direction(dir: str):
    raise ValueError(f"Unknown direction {dir}")


def main():
    pathway = Path(parser().parse_args().pathway)
    pathway_folder = intermediate_directory / pathway

    processed_directory.mkdir(exist_ok=True)

    # Create our output folder within the pathway directory
    out_folder = processed_directory / "pathways" / pathway
    out_folder.mkdir(exist_ok=True, parents=True)

    nodes_file = pathway_folder / "nodes.txt"
    nodes_df = pd.read_csv(nodes_file, sep="\t")

    # a dictionary mapping gene -> Uniprot accession ID
    gene_to_uniprot = pd.Series(nodes_df["uniprot"].values, index=nodes_df["NODE"]).to_dict()

    # nodes
    nodes_uniprot = nodes_df[["uniprot"]]
    nodes_uniprot.to_csv(out_folder / "gs_nodes.txt", sep="\t", index=False, header=False)

    # edges
    edges_df = pd.read_csv(pathway_folder / "edges.txt", sep="\t", header=0)
    edges_df = edges_df.rename(columns={"NODE1": "Interactor1", "NODE2": "Interactor2"})
    edges_df["Interactor1"] = edges_df["Interactor1"].map(gene_to_uniprot)
    edges_df["Interactor2"] = edges_df["Interactor2"].map(gene_to_uniprot)
    edges_df["Rank"] = 1
    edges_df["Direction"] = edges_df["INTERACTION_TYPE"].apply(
        lambda x: "D" if x in directed else ("U" if x in undirected else raise_unknown_direction(x))
    )
    edges_df = edges_df.drop(columns="INTERACTION_TYPE")

    # remove duplicate rows
    # sort by (node1 and node2) to ensure deterministic sorting
    edges_df = edges_df.sort_values(by=["Interactor1", "Interactor2"], ascending=True, ignore_index=True)
    undirected_mask = edges_df["Direction"] == "U"
    min_nodes = edges_df.loc[undirected_mask, ["Interactor1", "Interactor2"]].min(axis=1)
    max_nodes = edges_df.loc[undirected_mask, ["Interactor1", "Interactor2"]].max(axis=1)
    edges_df.loc[undirected_mask, "Interactor1"] = min_nodes
    edges_df.loc[undirected_mask, "Interactor2"] = max_nodes

    # keep 1 directed and 1 undirected edge if both exist
    # since rank is 1, we don't need to sort by rank.
    edges_df = edges_df.sort_values(by=["Interactor1", "Interactor2", "Direction"], ascending=True, ignore_index=True)
    edges_df = edges_df.drop_duplicates(keep="first", ignore_index=True)
    # We trim the gold standard edges against the interactome
    interactome_df = pd.read_csv(
        processed_directory / "interactome.tsv",
        sep="\t",
        header=None,
        names=["Interactor1", "Interactor2", "Weight", "Direction"],
        dtype={"Interactor1": str, "Interactor2": str},
    )
    edges_df = edges_df.merge(interactome_df, how="inner", on=["Interactor1", "Interactor2"])
    # We don't care about extraneous information provided by the interactome.
    edges_df = edges_df.drop(columns=["Direction_y", "Weight"]).rename(columns={"Direction_x": "Direction"})
    edges_df.to_csv(out_folder / "gs_edges.txt", sep="\t", index=False, header=False)

    # prizes, targets, sources
    prizes_file = pathway_folder / "prizes.txt"
    prizes_df = pd.read_csv(prizes_file, sep="\t")

    target_file = pathway_folder / "targets.txt"
    target_df = pd.read_csv(target_file, sep="\t")

    source_file = pathway_folder / "sources.txt"
    source_df = pd.read_csv(source_file, sep="\t")

    # final resulting df combining all the sources, targets, and prizes
    prizes_df["sources"] = prizes_df["uniprot"].isin(source_df["uniprot"])
    prizes_df["targets"] = prizes_df["uniprot"].isin(target_df["uniprot"])
    prizes_df["dummy"] = ""
    prizes_df.rename(columns={"uniprot": "NODEID", "prizes": "prize"}, inplace=True)
    result_df = prizes_df[["NODEID", "prize", "sources", "targets", "active", "dummy"]]
    # We trim the data file against the gold standard (which was already trimmed against the interactome)
    data_df = trim_data_file(data_df=result_df, gold_standard_df=edges_df)
    data_df.to_csv(out_folder / "node_prizes.txt", sep="\t", index=False, header=True)


if __name__ == "__main__":
    main()
