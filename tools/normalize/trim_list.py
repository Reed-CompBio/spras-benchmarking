import pandas
from pathlib import Path
from tools.normalize.interactome import get_interactome_nodes


def trim_node_list(interactome_nodes: list[str], nodes_to_trim: list[str]) -> list[str]:
    """
    Trims a node list with the desired interactome. The actual method you may want to use is in
    `trim_node_list_file`.
    """
    return list(set(nodes_to_trim).intersection(interactome_nodes))


def trim_node_list_file(interactome_df: pandas.DataFrame, node_list: Path, output: Path):
    """
    Trims a node list file with the desired interactome.0
    """
    return output.write_text(
        "\n".join(
            trim_node_list(
                list(get_interactome_nodes(interactome_df)), [line.strip() for line in node_list.read_text().splitlines() if line.strip() != ""]
            )
        )
    )
