import pandas


# def trim_input_nodes_file(interactome_df: pandas.DataFrame, input_nodes_df: pandas.DataFrame) -> pandas.DataFrame:
#     return input_nodes_df.loc[
#         input_nodes_df["NODEID"].isin(interactome_df["Interactor1"]) & input_nodes_df["NODEID"].isin(interactome_df["Interactor2"]), :
#     ]

def trim_input_nodes_file(interactome_df: pandas.DataFrame, input_nodes_df: pandas.DataFrame) -> pandas.DataFrame:
    """
    Keep input nodes that appear anywhere in the interactome (as Interactor1 or Interactor2)
    Nodes absent from the interactome cannot contribute to pathway reconstruction 
    """
    
    interactome_nodes = set(interactome_df["Interactor1"]).union(interactome_df["Interactor2"])
    return input_nodes_df.loc[input_nodes_df["NODEID"].isin(interactome_nodes), :]