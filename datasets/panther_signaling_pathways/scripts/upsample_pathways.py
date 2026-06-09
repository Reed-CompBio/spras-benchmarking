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
# - could we make more source focused pathway so  we add all pathways that are > 10 except for Cadherin%20signaling%20pathway and Wnt%20signaling%20pathway
# - could we make more balanced source and targets pathway and add pathways together to get better ratios of the source to target target to source
# - need to also pick which pathways to use as plain (I like to choose the Sources/Targets > 3 and Targets/Sources > 3 and see how many I am left with)
#   - but we can still use the other pathways for the upsampled ones

import sys
import pandas as pd
from pathlib import Path
panther_directory = Path(__file__).parent.parent.resolve()
intermediate_directory = panther_directory / "intermediate"
raw_directory = panther_directory / "raw"
processed_directory = panther_directory / "processed"


def main():
    
    stats = pd.read_csv(processed_directory / "filtered_stats.tsv", sep="\t")
    print(stats)
    pathways = stats["Name"].tolist()
    print(pathways)
    print(len(pathways))

    all_pathways = []
    all_input_nodes = []

    for p in pathways:
        print(p)
        pathway_folder = intermediate_directory / p
        pathway = pd.read_csv(pathway_folder / "pathway.csv", sep="\t")
        print(len(pathway))
        input_nodes = pd.read_csv(pathway_folder / "input_nodes.csv", sep="\t")
        all_pathways.append(pathway)
        all_input_nodes.append(input_nodes)

    # combine edges, drop duplicate edges
    mega_pathway = pd.concat(all_pathways, ignore_index=True)
    mega_pathway = mega_pathway.sort_values(
        by=["Node1", "Node2", "Direction"],
        ascending=[True, True, True],  # 'D' < 'U' so D stays first; highest weight first
    )
    mega_pathway = mega_pathway.drop_duplicates(subset=["Node1", "Node2"], keep="first").reset_index(drop=True)

    # combine input nodes, collapsing per NODEID
    # Stack every pathway's input-node rows into one frame; NODEID now repeats,
    # once per pathway that listed it.
    mega_input_nodes = pd.concat(all_input_nodes, ignore_index=True)

    # Per-column rule for collapsing those repeated NODEID rows back to one row.
    # "max" on a boolean acts as logical OR: a node flagged True in any pathway
    # stays True after the merge.
    agg = {
        "sources": "max",   # node is a source if it is a source in at least one pathway
        "targets": "max",   # node is a target if it is a target in at least one pathway
    }
    # prize and active are optional columns; include them in the merge rule only
    # if present. "max" keeps the largest prize seen for a node across pathways
    # and keeps active True if active in any pathway.
    if "prize" in mega_input_nodes.columns:
        agg["prize"] = "max"
    if "active" in mega_input_nodes.columns:
        agg["active"] = "max"

    # Collapse to one row per NODEID, applying the per-column rule above.
    # A node that is a source in one pathway and a target in another ends up
    # flagged as both.
    mega_input_nodes = mega_input_nodes.groupby("NODEID", as_index=False).agg(agg)

    mega_pathway.to_csv(processed_directory / "mega_pathway.csv", sep="\t", index=False)
    mega_input_nodes.to_csv(processed_directory / "mega_input_nodes.csv", sep="\t", index=False)
    print(f"Mega pathway: {len(mega_pathway)} edges, {len(mega_input_nodes)} input nodes")

        


if __name__ == "__main__":
    main()