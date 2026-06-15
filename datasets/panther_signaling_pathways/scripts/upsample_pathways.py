# TODO: Upsample not just all of the pathways together
# but some that are intermediate values between the pathways and the fully upsampled pathway
# TODO: Only change one at a time; upsample the pathway (and have intermediate upsampled pathways),
#  keep interactome the same (100%) and then downsample interactome and keep the pathways the same (100%)

# these are going to be used for the performance evaluation
# 1) upsample all the pathways together as one mega pathway
# - use filtered_stats.tsv to tell us what pathways to mega combine together that are useable
# - get the merged pathway and the merged input_nodes
# 2) want to have some intermediate points that are the values between for the sources and targets
# - individual source/target counts sum to roughly 350 source-slots and 175 target-slots,
# but the mega has only 327 total input nodes, so there's heavy sharing between the pathways

# - could we make a pathway that focuses more on targets so we add all the pathways together that have targets > 10 other than for Wnt%20signaling%20pathway
#   - one without the 1 and one with the 1
# - could we make more source focused pathway so  we add all pathways that are > 10 except for Cadherin%20signaling%20pathway and Wnt%20signaling%20pathway
#   - one without those 2 and one with those 2

# - could we make more balanced source and targets pathway and add pathways together to get better ratios of the source to target target to source
# - need to also pick which pathways to use as plain (I like to choose the Sources/Targets > 3 and Targets/Sources > 3 and see how many I am left with)
#   - but we can still use the other pathways for the upsampled ones


"""
upsample_pathways.py

Combines individual PANTHER signaling pathways into merged "upsampled" pathway
variants for benchmarking. Each variant is a different subset of the filtered
pathways, combined into a single edge list and input-node list.

Output files are written to processed/upsampled/:
    mega_pathway.csv / mega_input_nodes.csv  -- all pathways merged
    target_focused_pathway.csv / _input_nodes.csv  -- pathways with > TARGET_MIN targets, excluding Wnt
    source_focused_pathway.csv / ...     -- pathways with > SOURCE_MIN sources, excluding Cadherin and Wnt
    balanced_pathway.csv / ...    -- greedy-balanced subset (sources ~ targets)
    upsampled_stats.tsv     -- source/target counts for each variant
"""

import pandas as pd
from pathlib import Path

panther_directory   = Path(__file__).parent.parent.resolve()
intermediate_directory = panther_directory / "intermediate"
processed_directory    = panther_directory / "processed"
upsampled_directory    = processed_directory / "upsampled"


# URL-encoded names matching the folder names under intermediate/
CADHERIN = "Cadherin%20signaling%20pathway"
WNT = "Wnt%20signaling%20pathway"

# Thresholds for "source-rich" and "target-rich" pathway selection
TARGET_MIN = 10
SOURCE_MIN = 10


def load_pathway(name):
    """
    Load a pathway's edge list and input-node list from intermediate/
    """
    folder = intermediate_directory / name
    edges = pd.read_csv(folder / "pathway.csv", sep="\t")
    input_nodes = pd.read_csv(folder / "input_nodes.csv", sep="\t")
    return edges, input_nodes

def count_roles(input_nodes):
    """
    Returns amount of sources and targets
    """
    n_sources = int((input_nodes["sources"]).sum())
    n_targets = int((input_nodes["targets"]).sum())
    return n_sources, n_targets

def merge_edges(edge_frames):
    """
    Concatenate edge frames and deduplicate by (Node1, Node2).
    When two pathways have the same edge with different directions,
    the directed edge is kept over the undirected.
    """
    merged = pd.concat(edge_frames, ignore_index=True)
    merged = merged.sort_values(["Node1", "Node2", "Direction"], ascending=True)
    merged = merged.drop_duplicates(subset=["Node1", "Node2"], keep="first")
    return merged.reset_index(drop=True)


def merge_input_nodes(node_frames):
    """
    Concatenate input-node frames and collapse to one row per NODEID.
    A node is a source/target/active if it is a source/target/active in any pathway (logical OR).
    The prize column takes the max across pathways.
    """
    merged = pd.concat(node_frames, ignore_index=True)

    # "any" is logical OR; a node flagged True in any pathway stays True.
    # For numeric columns take the max.
    agg = {
        "sources": "any",
        "targets": "any",
        "prize":   "max",
        "active":  "any",
    }

    return merged.groupby("NODEID", as_index=False).agg(agg)

def merge_pathways(names, edges, nodes):
    """
    Merge edges and input nodes for the given list of pathway names.
    """
    merged_edges = merge_edges([edges[n] for n in names])
    merged_nodes = merge_input_nodes([nodes[n] for n in names])
    return merged_edges, merged_nodes

def write_variant(tag, names, edges, nodes):
    """
    Merge the pathways in names, write the results to upsampled/, and print a
    summary. Returns the merged input-node DataFrame (used for stats), or None
    if the name list is empty.
    """

    pathway_df, input_nodes_df = merge_pathways(names, edges, nodes)
    n_sources, n_targets = count_roles(input_nodes_df)

    pathway_df.to_csv(upsampled_directory / f"{tag}_pathway.csv",     sep="\t", index=False)
    input_nodes_df.to_csv(upsampled_directory / f"{tag}_input_nodes.csv", sep="\t", index=False)

    print(f"{tag}: {len(names)} pathways -> {len(pathway_df)} edges, "
          f"{len(input_nodes_df)} nodes ({n_sources} sources, {n_targets} targets)")
    return input_nodes_df


def pathway_stats_row(name, input_nodes):
    """
    Build one row for the summary stats table.
    """
    s, t = count_roles(input_nodes)
    return {
        "Name": name,
        "Sources": s,
        "Targets": t,
        "Ratio Sources/Targets": s / t if t != 0 else float("nan"),
        "Ratio Targets/Sources": t / s if s != 0 else float("nan"),
    }

def imbalance(n_sources, n_targets):
    """
    Symmetric imbalance score.
    Returns max/min - 1

    If sources=100 and targets=125, then 125/100 - 1 = 0.25,
    meaning the larger side is 25% bigger than the smaller.
    A tolerance of 0.25 allows up to that ratio, and 0.0 would require exact equality.

    """
    return max(n_sources, n_targets) / min(n_sources, n_targets) - 1.0

# Balanced-subset selection
def grow_balanced(pathway_names, edges, nodes, tol=0.25):
    """
    Build a balanced subset by greedy expansion.

    1. Seed with the single most balanced pathway.
    2. At each step, try adding each remaining pathway. Pick the one that
       maximizes merged input-node count while keeping imbalance within `tol`.
    3. Stop when no remaining pathway keeps the merge within tolerance.

    Returns selected_names.
    """
    pool = list(pathway_names)

    # Seed: the single pathway closest to a 1:1 source/target ratio
    seed = min(pool, key=lambda n: imbalance(*count_roles(nodes[n])))
    selected = [seed]
    pool.remove(seed)

    # TODO: what if we randomize the pool? right now this is based on the order of the pool
    while pool:
        best_candidate = None
        best_node_count = len(merge_input_nodes([nodes[n] for n in selected]))

        for candidate in pool:
            merged_nodes = merge_input_nodes([nodes[n] for n in selected + [candidate]])
            cs, ct = count_roles(merged_nodes)

            if imbalance(cs, ct) <= tol and len(merged_nodes) > best_node_count:
                best_candidate  = candidate
                best_node_count = len(merged_nodes)

        if best_candidate is None:
            break  # no candidate improves size while staying balanced

        selected.append(best_candidate)
        pool.remove(best_candidate)

    return selected

def main():
    upsampled_directory.mkdir(parents=True, exist_ok=True)

    # Load the list of filtered pathways
    stats = pd.read_csv(processed_directory / "filtered_stats.tsv", sep="\t")
    pathway_names = stats["Name"].tolist()

    # Load all pathways once; reuse across all variants
    edges, nodes, role_counts = {}, {}, {}
    for name in pathway_names:
        edges[name], nodes[name] = load_pathway(name)
        role_counts[name]        = count_roles(nodes[name])

    stats_rows = []

    #  Mega: all pathways merged
    print("\nBuilding variants:")
    mega_edges, mega_nodes = merge_pathways(pathway_names, edges, nodes)
    mega_edges.to_csv(upsampled_directory / "mega_pathway.csv",     sep="\t", index=False)
    mega_nodes.to_csv(upsampled_directory / "mega_input_nodes.csv", sep="\t", index=False)
    s, t = count_roles(mega_nodes)
    print(f"mega: {len(pathway_names)} pathways -> {len(mega_edges)} edges, "
          f"{len(mega_nodes)} nodes ({s} sources, {t} targets)")
    stats_rows.append(pathway_stats_row("mega", mega_nodes))

    # Source or Target variants
    target_rich = [n for n in pathway_names if role_counts[n][1] > TARGET_MIN]
    source_rich = [n for n in pathway_names if role_counts[n][0] > SOURCE_MIN]

    variants = [
        # Target-focused: pathways with > TARGET_MIN targets
        # Wnt clears the cutoff but floods the source side, so we remove it
        ("target_focused", [n for n in target_rich if n != WNT]),

        # Source-focused: pathways with > SOURCE_MIN sources
        # Cadherin and Wnt both dominate the source side, so we remove them
        ("source_focused",[n for n in source_rich if n not in (CADHERIN, WNT)]),
    ]

    for tag, names in variants:
        variant_nodes = write_variant(tag, names, edges, nodes)
        if variant_nodes is not None:
            stats_rows.append(pathway_stats_row(tag, variant_nodes))

    # Balanced variant
    # tol=0.15 is stricter;
    for tol, tag in [(0.15, "balanced")]:
        balanced_names = grow_balanced(pathway_names, edges, nodes, tol=tol)
        print(f"{tag} (tol {tol})")
        print(f"pathways: {balanced_names}")
        balanced_nodes = write_variant(tag, balanced_names, edges, nodes)
        if balanced_nodes is not None:
            stats_rows.append(pathway_stats_row(tag, balanced_nodes))

    # Write summary stats
    pd.DataFrame(stats_rows).to_csv(upsampled_directory / "upsampled_stats.tsv", sep="\t", index=False)
    print(f"\nWrote stats for {len(stats_rows)} variants to upsampled_stats.tsv")


if __name__ == "__main__":
    main()