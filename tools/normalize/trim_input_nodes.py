import pandas


def trim_input_nodes_file(interactome_df: pandas.DataFrame, input_nodes_df: pandas.DataFrame) -> pandas.DataFrame:
    return input_nodes_df.loc[
        input_nodes_df["NODEID"].isin(interactome_df["Interactor1"]) & input_nodes_df["NODEID"].isin(interactome_df["Interactor2"]), :
    ]
