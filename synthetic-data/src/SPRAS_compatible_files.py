import os
import pandas as pd

spras_compatible_dir = "spras-compatible-pathway-data/"
if not os.path.exists(spras_compatible_dir):
    os.makedirs(spras_compatible_dir)

pathway_dirs = ["Apoptosis_signaling", "B_cell_activation", "Beta3_adrenergic_rec", "Cadherin_signaling", "Hedgehog_signaling", "Insulin_IGF", "Interleukin_signaling", "Notch_signaling", "PDGF_signaling", "Ras", "T_cell_activation", "Toll_signaling", "Wnt_signaling", "p38_MAPK", "Nicotinic_acetylchol"]
directory = "pathway-data/"


for pathway in pathway_dirs:
    pathway_folder = directory + pathway + "/"

    # Create the output folder "uniprot" within the pathway directory
    out_folder = os.path.join(spras_compatible_dir, pathway)
    print(out_folder)
    os.makedirs(out_folder, exist_ok=True)
    
    nodes_file = os.path.join(pathway_folder, "NODES.txt")
    nodes_df = pd.read_csv(nodes_file, sep="\t")
    
    # a dictionary mapping gene -> Uniprot accession ID
    gene_to_uniprot = pd.Series(nodes_df['uniprot'].values, index=nodes_df['NODE']).to_dict()
    
    nodes_uniprot = nodes_df[['uniprot']]
    nodes_uniprot.to_csv(os.path.join(out_folder, f"{pathway}_gs_nodes.txt"), sep="\t", index=False, header=False)

    edges_file = os.path.join(pathway_folder, "EDGES.txt")
    edges_df = pd.read_csv(edges_file, sep="\t")
    edges_df['NODE1'] = edges_df['NODE1'].map(gene_to_uniprot)
    edges_df['NODE2'] = edges_df['NODE2'].map(gene_to_uniprot)
    edges_df['Rank'] = 1
    edges_df['Direction'] = "U"
    edges_df.to_csv(os.path.join(out_folder, f"{pathway}_gs_edges.txt"), sep="\t", index=False, header=False)

    prizes_file = os.path.join(pathway_folder, "PRIZES-100.txt")
    prizes_df = pd.read_csv(prizes_file, sep="\t")
    prizes_uniprot = prizes_df[['uniprot', 'prizes', 'active']]

    target_file = os.path.join(pathway_folder, "TARGETS.txt")
    target_df = pd.read_csv(target_file, sep="\t")
    target_uniprot = target_df[['uniprot']]
    
    source_file = os.path.join(pathway_folder, "SOURCES.txt")
    source_df = pd.read_csv(source_file, sep="\t")
    source_uniprot = source_df[['uniprot']]
    
    # final resulting df combining all the sources, targets, and prizes
    prizes_df['sources'] = prizes_df['uniprot'].isin(source_df['uniprot'])
    prizes_df['targets'] = prizes_df['uniprot'].isin(target_df['uniprot'])
    prizes_df['dummy'] = ""
    prizes_df.rename(columns={'uniprot': 'NODEID', 'prizes': 'prize'}, inplace=True)
    result_df = prizes_df[['NODEID', 'prize', 'sources', 'targets', 'active', 'dummy']]
    result_df.to_csv(os.path.join(out_folder, f"{pathway}_node_prizes.txt"), sep="\t", index=False, header=True)
