import pandas
from pathlib import Path

current_directory = Path(__file__).parent.resolve()

def main():
    # convert the interactome to SPRAS format
    print('Reading interactome...')
    interactome_df = pandas.read_csv(current_directory / '..' / 'raw' / '9606.protein.links.full.v12.0.txt', sep=' ',
                         usecols=['protein1', 'protein2', 'combined_score'])
    interactome_df.columns = ['Protein1', 'Protein2', 'Weight']

    # we also want to representatively remove a certain percentage of elements from the interactome,
    # to make sure our interactome downsampling preserves edge weight distributions
    # (we don't care to preserve other major topological properties just yet.)
    # since this file is large, we opt for streaming the interactome for removing edges instead
    
    print('Initially processing interactome...')
    interactome_df['Weight'] = interactome_df['Weight'].div(1000) # scores are from 1-1000: we normalize from 0-1.
    interactome_df['Direction'] = 'U'
    print('Sorting interactome...')
    interactome_df = interactome_df.sort_values('Weight', kind='stable')

    print('Mapping interactome...')
    # STRINGDB -> UniProt accession ID pairings
    UniProt_AC = pandas.read_csv(current_directory / '..' / "raw" / "human-interactome" / "String_to_Uniprot_ids_2025_04_06.tsv", sep='\t', header=0)
    one_to_many_dict = UniProt_AC.groupby("From")["Entry"].apply(list).to_dict()
    def get_aliases(protein_id):
        return one_to_many_dict.get(protein_id, [])

    interactome_df['Protein1_uniprot'] = interactome_df['Protein1'].apply(get_aliases)
    interactome_df['Protein2_uniprot'] = interactome_df['Protein2'].apply(get_aliases)

    interactome_df = interactome_df.explode('Protein1_uniprot').explode('Protein2_uniprot')

    missing_alias_edges = interactome_df[(interactome_df['Protein1_uniprot'].isna()) | (interactome_df['Protein2_uniprot'].isna())]

    proteins_without_aliases = pandas.concat([
        missing_alias_edges.loc[missing_alias_edges['Protein1_uniprot'].isna(), 'Protein1'],
        missing_alias_edges.loc[missing_alias_edges['Protein2_uniprot'].isna(), 'Protein2']
    ], ignore_index=True).drop_duplicates().reset_index(drop=True)
    proteins_without_aliases = proteins_without_aliases.to_frame(name="protein")

    removed_edges = missing_alias_edges[['Protein1', 'Protein2']]
    removed_edges = removed_edges.drop_duplicates().reset_index(drop=True)

    (current_directory / '..' / "processed" / "interactomes" / "uniprot-threshold-interactomes").mkdir(exist_ok=True, parents=True)
    proteins_without_aliases.to_csv(
        current_directory / '..' / "processed" / "interactomes" / "uniprot-threshold-interactomes" / "proteins_missing_aliases.csv",
        sep='\t', index=False, header=True)
    removed_edges.to_csv(
        current_directory / '..' / "processed" / "interactomes" / "uniprot-threshold-interactomes" / "removed_edges.txt",
        sep='\t', index=False, header=True)
    interactome_df = interactome_df.dropna(subset=['Protein1_uniprot', 'Protein2_uniprot']).reset_index(drop=True)
    interactome_df = interactome_df[["Protein1_uniprot", "Protein2_uniprot", "Weight", "Direction"]]

    print('Counting weight counts...')
    interactome_df['Weight'].value_counts(sort=False).to_csv(
        current_directory / '..' / 'processed' / 'weight-counts.tsv', sep='\t')

    print('Saving interactome...')
    interactome_df.to_csv(current_directory / '..' / 'processed' / 'interactome.tsv', sep='\t', header=False, index=False)


if __name__ == "__main__":
    main()
