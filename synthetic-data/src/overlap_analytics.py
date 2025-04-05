import pandas as pd
import glob

human_interactome = pd.read_csv("preprocessing/human-interactome/9606.protein.links.full.v12.0.txt", sep = " ", header=0)
human_interactome.sort_values("experiments_transferred", inplace=True)

thresholds = [1, 100, 200, 300, 400, 500, 600, 700, 800, 900]
pathway_dirs = ["Apoptosis_signaling", "B_cell_activation", "Beta3_adrenergic_rec", "Cadherin_signaling", "Hedgehog_signaling", "Insulin_IGF", "Interleukin_signaling", "Notch_signaling", "PDGF_signaling", "Ras", "T_cell_activation", "Toll_signaling", "Wnt_signaling", "p38_MAPK"]

# overlap of edges in pathway per threshold
for threshold in thresholds:
    threshold_human_interactome = human_interactome[
        (human_interactome["experiments"] >= threshold) | 
        (human_interactome["experiments_transferred"] >= threshold)
    ]

    print(f"threshold {threshold}; # of edges {len(threshold_human_interactome)}")

    for pathway_dir in pathway_dirs:
        files = glob.glob(f"preprocessing/{pathway_dir}/STRING-EDGES*")[0]
        edges = pd.read_csv(files, sep = "\t", header=0)
        print(f"{pathway_dir}; # of edges = {len(edges)}")
        edges.columns = ["protein1", "protein2"]
        
        overlap = pd.merge(threshold_human_interactome, edges, on=['protein1', 'protein2'])
        print(f"# of common edges between {pathway_dir} pathway and interactome: {len(overlap)}")
        percent = len(overlap) / len(edges)
        print(f"% of edges included = {percent}")

    print()
    print()

# overall overlap of all the edges
combined_edges = pd.DataFrame(columns=["protein1", "protein2"])
for pathway_dir in pathway_dirs:
    files = glob.glob(f"preprocessing/{pathway_dir}/STRING-EDGES*")[0]
    edges = pd.read_csv(files, sep = "\t", header=0)
    edges.columns = ["protein1", "protein2"]
    combined_edges = pd.concat([combined_edges, edges], ignore_index=True)

combined_edges.drop_duplicates(inplace=True)
print(f"Total combined edges: {len(combined_edges)}")

for threshold in thresholds:
    threshold_human_interactome = human_interactome[
        (human_interactome["experiments"] >= threshold) | 
        (human_interactome["experiments_transferred"] >= threshold)
    ]

    print(f"threshold {threshold}; # of edges {len(threshold_human_interactome)}")

    overlap = pd.merge(threshold_human_interactome, combined_edges, on=['protein1', 'protein2'])
    print(f"# of common edges between pathways and interactome: {len(overlap)}")
    percent = len(overlap) / len(combined_edges)
    print(f"% of edges included from pathway = {percent}")


# put in a df to be more readable
results = []

for threshold in thresholds:
    threshold_human_interactome = human_interactome[
        (human_interactome["experiments"] >= threshold) | 
        (human_interactome["experiments_transferred"] >= threshold)
    ]
    interactome_count = len(threshold_human_interactome)
    
    for pathway_dir in pathway_dirs:
        files = glob.glob(f"preprocessing/{pathway_dir}/STRING-EDGES*")[0]
        edges = pd.read_csv(files, sep="\t", header=0)
        edges.columns = ["protein1", "protein2"]
        edges_count = len(edges)
        
        overlap = pd.merge(threshold_human_interactome, edges, on=['protein1', 'protein2'])
        overlap_count = len(overlap)
        percent_overlap = overlap_count / edges_count if edges_count != 0 else 0

        results.append({
            "threshold": threshold,
            "pathway": pathway_dir,
            "#_of_interactome_edges": interactome_count,
            "pathway_edges_count": edges_count,
            "overlap_between_pathway_and_interactome_count": overlap_count,
            "percentage_of_edges_included": percent_overlap
        })

overlap_info_df = pd.DataFrame(results)

overlap_info_df.to_csv("overlap_info_df.csv", sep="\t", index=False)