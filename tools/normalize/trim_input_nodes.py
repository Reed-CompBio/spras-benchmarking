import pandas


def trim_input_nodes_file(interactome_df: pandas.DataFrame, input_nodes_df: pandas.DataFrame) -> pandas.DataFrame:
    """
    Keep input nodes that appear anywhere in the interactome.

    Nodes absent from the interactome cannot contribute to pathway reconstruction.

    Both inputs must be in SPRAS format. The interactome is expected to have the
    columns Interactor1, Interactor2, Weight, Direction.
    """
    interactome_nodes = set(interactome_df["Interactor1"]).union(interactome_df["Interactor2"])
    return input_nodes_df.loc[input_nodes_df["NODEID"].isin(interactome_nodes), :]