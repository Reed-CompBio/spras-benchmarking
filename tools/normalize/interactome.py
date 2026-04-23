import pandas

def deduplicate_edges(interactome_df: pandas.DataFrame) -> tuple[pandas.DataFrame, bool]:
    """
    Removes duplicate edges from the input `DataFrame` as an interactome.
    - For duplicate edges (based on Node1, Node2, and Direction), the one with the smallest Rank is kept.
    - For undirected edges, the node pair is sorted (e.g., "B-A" becomes "A-B") before removing duplicates.

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

    unique_edges_df = interactome_df.drop_duplicates(subset=["Interactor1", "Interactor2", "Direction"], keep="first", ignore_index=True)

    return unique_edges_df, not unique_edges_df.equals(interactome_df)

def get_interactome_nodes(interactome_df: pandas.DataFrame) -> set[str]:
    """
    Gets all nodes associated with an interactome.

    @param interactome_df: The pandas dataframe with Interactor1 and Interactor2 as columns.
    @return the list of nodes in Interactor1 and Interactor2.
    
    NOTE: This isn't guaranteed to be order stable, not for any externally-required reason.
    """
    return set(list(interactome_df["Interactor1"])).union(list(interactome_df["Interactor2"]))

# TODO: implement in SPRAS? This is the direct analogue of duplicate_edges in util.
# We should also separate deduplication from the broader normalization if we do any more steps.
def normalize_interactome(interactome_df: pandas.DataFrame) -> tuple[pandas.DataFrame, bool]:
    """
    Applies a number of post-processing steps to an interactome. This currently only includes:
    - Deduplicating undirected edges
    """
    return deduplicate_edges(interactome_df)
