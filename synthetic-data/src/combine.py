import pandas as pd
import os

if not os.path.exists("interactomes/"):
    os.makedirs("interactomes/")

if not os.path.exists("interactomes/uniprot-combined-threshold-interactomes/"):
    os.makedirs("interactomes/uniprot-combined-threshold-interactomes/")

thresholds = [1, 100, 200, 300, 400, 500, 600, 700, 800, 900]
pathway_dirs = ["Apoptosis_signaling", "B_cell_activation", "Beta3_adrenergic_rec", "Cadherin_signaling", "Hedgehog_signaling", "Insulin_IGF", "Interleukin_signaling", "Notch_signaling", "PDGF_signaling", "Ras", "T_cell_activation", "Toll_signaling", "Wnt_signaling", "p38_MAPK", "Nicotinic_acetylchol"]

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