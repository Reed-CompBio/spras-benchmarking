import pandas

def trim_data_file(data_df: pandas.DataFrame, gold_standard_df: pandas.DataFrame) -> pandas.DataFrame:
    """
    Trims the associated SPRAS data file with the nodes in the gold standard file.
    """
    # We just want the set of all nodes present in the gold standard
    gold_standard_nodes = set(gold_standard_df["Interactor1"]).union(set(gold_standard_df["Interactor2"]))
    return data_df[data_df["NODEID"].isin(gold_standard_nodes)]
