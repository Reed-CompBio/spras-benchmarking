import pandas
from pathlib import Path

current_directory = Path(__file__).parent.resolve()

interactome_folder = current_directory / ".." / "raw" / "human-interactome"

def main():
    tf_df = pandas.read_csv(interactome_folder / "Homo_sapiens_TF.txt", sep = "\t", header = 0)
    df = pandas.read_csv(interactome_folder / 'Ensembl_to_Uniprot_ids_2025_04_08.tsv', sep='\t', header = 0)

    def filter_group(group):
        if 'reviewed' in group['Reviewed'].values:
            return group[group['Reviewed'] == 'reviewed']
        else:
            return group

    filtered_df = df.groupby(by = ['From'], group_keys=False).apply(filter_group)
    filtered_df.to_csv(interactome_folder / "filtered_Ensembl_to_Uniprot_ids_2025_04_08.txt", sep = '\t', header = True, index = False)

    one_to_many_dict = filtered_df.groupby("From")["Entry"].apply(list).to_dict()

    def get_aliases(protein_id):
        return one_to_many_dict.get(protein_id, [])

    tf_df['Uniprot_Accession'] = tf_df['Ensembl'].apply(get_aliases)
    tf_df = tf_df.explode('Uniprot_Accession')
    tf_df = tf_df.fillna('NA')
    tf_df.to_csv(interactome_folder / "Homo_sapiens_TF_Uniprot.txt", header = True, sep = "\t", index = False)

if __name__ == "__main__":
    main()
