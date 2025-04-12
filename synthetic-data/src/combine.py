import pandas as pd
import os

if not os.path.exists("interactomes/"):
    os.makedirs("interactomes/")

if not os.path.exists("interactomes/uniprot-combined-threshold-interactomes/"):
    os.makedirs("interactomes/uniprot-combined-threshold-interactomes/")

thresholds = [1, 100, 200, 300, 400, 500, 600, 700, 800, 900]
pathway_dirs = ["Apoptosis_signaling", "B_cell_activation", "Beta3_adrenergic_rec", "Cadherin_signaling", "Hedgehog_signaling", "Insulin_IGF", "Interleukin_signaling", "Notch_signaling", "PDGF_signaling", "Ras", "T_cell_activation", "Toll_signaling", "Wnt_signaling", "p38_MAPK", "Nicotinic_acetylchol", "Fas_signaling", "FGF_signaling", "Interferon_gamma_signaling", "JAK_STAT_signaling", "VEGF_signaling"] 

# pilot data
# ["Wnt_signaling", "JAK_STAT_signaling", "Interferon_gamma_signaling", "FGF_signaling", "Ras" ]

combined_edges = pd.DataFrame(columns=["Node1", "Node2"])
for pathway_dir in pathway_dirs:
    edge_file = f"spras-compatible-pathway-data/{pathway_dir}/{pathway_dir}_gs_edges.txt"
    edges = pd.read_csv(edge_file, sep = "\t")
    edges.columns = ["Node1", "Node2", "Rank", "Direction"]
    combined_edges = pd.concat([combined_edges, edges], ignore_index=True)

combined_edges.drop_duplicates(inplace=True)

for threshold in thresholds:
    threshold_human_interactome = pd.read_csv(f"interactomes/uniprot-threshold-interactomes/uniprot_human_interactome_{threshold}.txt", sep="\t")
    threshold_human_interactome.columns = ["Node1", "Node2", "Rank", "Direction"]

    fifty_percentile_rank = threshold_human_interactome.describe().loc['50%', "Rank"]
    combined_edges["Rank"] = fifty_percentile_rank

    merged_df = pd.concat([combined_edges, threshold_human_interactome])
    merged_df = merged_df.drop_duplicates(subset=['Node1', 'Node2'], keep='first')
    
    merged_df.to_csv(f"interactomes/uniprot-combined-threshold-interactomes/uniprot_combined_interactome_{threshold}.txt", sep="\t", index = False, header = False)

# get overlap analytics

# overlap of edges in pathway per threshold
results = []

for threshold in thresholds:
    threshold_human_interactome = pd.read_csv(f"interactomes/uniprot-threshold-interactomes/uniprot_human_interactome_{threshold}.txt", sep="\t")
    threshold_human_interactome.columns = ["Node1", "Node2", "Rank", "Direction"]

    interactome_count = len(threshold_human_interactome)
    
    for pathway_dir in pathway_dirs:
        edge_file = f"spras-compatible-pathway-data/{pathway_dir}/{pathway_dir}_gs_edges.txt"
        edges = pd.read_csv(edge_file, sep = "\t")
        edges.columns = ["Node1", "Node2", "Rank", "Direction"]
        edges_count = len(edges)
        
        overlap = pd.merge(threshold_human_interactome, edges, on=['Node1', 'Node2'])
        overlap_count = len(overlap)
        percent_overlap = overlap_count / edges_count if edges_count != 0 else 0

        edges_to_add = edges_count - overlap_count

        results.append({
            "threshold": threshold,
            "pathway": pathway_dir,
            "number_of_interactome_edges": interactome_count,
            "number_of_pathway_edges": edges_count,
            "overlap_between_pathway_and_interactome_count": overlap_count,
            "percentage_of_edges_included": percent_overlap,
            "number_of_edges_to_add": edges_to_add
        })

overlap_info_df = pd.DataFrame(results)
overlap_info_df.to_csv("interactomes/uniprot-combined-threshold-interactomes/overlap_info.csv", sep="\t", index=False)

# overall overlap of all the edges
results = []

combined_edges = pd.DataFrame(columns=["Node1", "Node2"])
for pathway_dir in pathway_dirs:
    edge_file = f"spras-compatible-pathway-data/{pathway_dir}/{pathway_dir}_gs_edges.txt"
    edges = pd.read_csv(edge_file, sep = "\t")
    edges.columns = ["Node1", "Node2", "Rank", "Direction"]
    combined_edges = pd.concat([combined_edges, edges], ignore_index=True)
combined_edges.drop_duplicates(inplace=True)
total_combined_edges = len(combined_edges)

for threshold in thresholds:
    interactome_file = f"interactomes/uniprot-threshold-interactomes/uniprot_human_interactome_{threshold}.txt"
    threshold_human_interactome = pd.read_csv(interactome_file, sep="\t")
    threshold_human_interactome.columns = ["Node1", "Node2", "Rank", "Direction"]
    interactome_count = len(threshold_human_interactome)

    # Find overlapping edges between the interactome and combined pathway edges
    overlap = pd.merge(threshold_human_interactome, combined_edges, on=['Node1', 'Node2'])
    overlap_count = len(overlap)
    
    # Calculate the percentage of pathway edges included in the interactome overlap
    percent_overlap = overlap_count / total_combined_edges if total_combined_edges > 0 else 0
    
    edges_to_add = total_combined_edges - overlap_count

    # Append the information as a dictionary to the results list
    results.append({
        "threshold": threshold,
        "total_interactome_edges": interactome_count,
        "total_combined_pathway_edges": total_combined_edges,
        "overlap_count": overlap_count,
        "percentage_of_edges_included": percent_overlap,
        "number_of_edges_to_add": edges_to_add
    })

# Create a DataFrame from the results list
overlap_combined_info_df = pd.DataFrame(results)
overlap_combined_info_df.to_csv("interactomes/uniprot-combined-threshold-interactomes/overlap_combined_info.csv", sep="\t", index=False)
