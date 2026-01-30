import pandas as pd

tf = pd.read_csv("human-interactome/Homo_sapiens_TF.txt", sep = "\t", header = 0)

df = pd.read_csv('human-interactome/Ensembl_to_Uniprot_ids_2025_04_08.tsv', sep='\t', header = 0)

def filter_group(group):
    if 'reviewed' in group['Reviewed'].values:
        return group[group['Reviewed'] == 'reviewed']
    else:
        return group

filtered_df = df.groupby(by = ['From'], group_keys=False).apply(filter_group)
filtered_df.to_csv("human-interactome/filtered_Ensembl_to_Uniprot_ids_2025_04_08.txt", sep = '\t', header = True, index = False)

one_to_many_dict = filtered_df.groupby("From")["Entry"].apply(list).to_dict()

def get_aliases(protein_id):
    return one_to_many_dict.get(protein_id, [])

tf['Uniprot_Accession'] = tf['Ensembl'].apply(get_aliases)
tf = tf.explode('Uniprot_Accession')
tf = tf.fillna('NA')
tf.to_csv("human-interactome/Homo_sapiens_TF_Uniprot.txt", header = True, sep = "\t", index = False)