import pandas as pd


def deduplicate_edges(interactome_df: pd.DataFrame) -> tuple[pd.DataFrame, bool]:
    """
    Removes duplicate edges from the input `DataFrame` as an interactome.
    - For undirected edges, the node pair is sorted (e.g., "B-A" becomes "A-B") before removing duplicates.
    - Directed edges are deduplicated by exact (Interactor1, Interactor2) match; "A->B" and "B->A" remain distinct.
    - When duplicates are found, the row with the highest Weight is retained.

    @param interactome_df: A DataFrame from a raw file pathway.
    @return pd.DataFrame: A DataFrame with duplicate edges removed.
    @return bool: True if duplicate edges were found and removed, False otherwise.
    """
    # for undirected edges, sort node pairs so that Node1 is always the lesser of the two
    undirected_mask = interactome_df["Direction"] == "U"

    # computes the minimum and maximum node (sorted order) for each row under the mask
    min_nodes = interactome_df.loc[undirected_mask, ["Interactor1", "Interactor2"]].min(axis=1)
    max_nodes = interactome_df.loc[undirected_mask, ["Interactor1", "Interactor2"]].max(axis=1)

    # assigns the sorted Node1 and Node2 back into the df
    interactome_df.loc[undirected_mask, "Interactor1"] = min_nodes
    interactome_df.loc[undirected_mask, "Interactor2"] = max_nodes

    # sort by Weight descending so keep="first" retains the highest-weight duplicate
    interactome_df.sort_values("Weight", ascending=False, inplace=True)

    # remove duplicated edges
    unique_edges_df = interactome_df.drop_duplicates(subset=["Interactor1", "Interactor2", "Direction"], keep="first", ignore_index=True)

    # return True if duplicates edges are removed
    duplicates_removed = len(unique_edges_df) < len(interactome_df)

    return unique_edges_df, duplicates_removed


def get_interactome_nodes(interactome_df: pd.DataFrame) -> set[str]:
    """
    Gets all nodes associated with an interactome.

    @param interactome_df: The pd dataframe with Interactor1 and Interactor2 as columns.
    @return the list of nodes in Interactor1 and Interactor2.

    NOTE: This isn't guaranteed to be order stable, not for any externally-required reason.
    """
    return set(list(interactome_df["Interactor1"])).union(list(interactome_df["Interactor2"]))
