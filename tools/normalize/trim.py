import pandas as pd
from pathlib import Path
from tools.normalize.interactome import get_interactome_nodes


def trim_node_list(interactome_nodes: list[str], nodes_to_trim: list[str]) -> list[str]:
    """
    Trims a node list with the desired interactome. The actual method you may want to use is in
    `trim_node_list_file`.
    """
    return list(set(nodes_to_trim).intersection(interactome_nodes))


def trim_node_list_file(interactome_df: pd.DataFrame, node_list: Path, output: Path):
    """
    Trims a node list file with the desired interactome.
    """
    return output.write_text(
        "\n".join(
            trim_node_list(
                list(get_interactome_nodes(interactome_df)), [line.strip() for line in node_list.read_text().splitlines() if line.strip() != ""]
            )
        )
    )

def trim_input_nodes_file(interactome_df: pd.DataFrame, input_nodes_df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep input nodes that appear anywhere in the interactome.

    Nodes absent from the interactome cannot contribute to pathway reconstruction.

    Both inputs must be in SPRAS format. The interactome is expected to have the
    columns Interactor1, Interactor2, Weight, Direction as the header. The input nodes is expected to at least have one column of NODEID as the header.
    """
    interactome_nodes = set(interactome_df["Interactor1"]).union(interactome_df["Interactor2"])
    return input_nodes_df.loc[input_nodes_df["NODEID"].isin(interactome_nodes), :]

def trim_gold_standard_nodes_file(interactome_df: pd.DataFrame, gs_df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep gold standard nodes that appear anywhere in the interactome.

    Nodes absent from the interactome cannot contribute to pathway reconstruction.

    Both inputs must be in SPRAS format. The interactome is expected to have the
    columns Interactor1, Interactor2, Weight, Direction as the header. The gold standard nodes is expected to have one column of NODEID as the header.
    """
    interactome_nodes = set(interactome_df["Interactor1"]).union(interactome_df["Interactor2"])
    return gs_df.loc[gs_df["NODEID"].isin(interactome_nodes), :]