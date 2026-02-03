import pandas
from pathlib import Path

current_directory = Path(__file__).parent.resolve()

interactome_folder = current_directory / ".." / "raw" / "human-interactome"

def main():
    tf_df = pandas.read_csv(interactome_folder / "Homo_sapiens_TF.tsv", sep = "\t", header = 0)
    # The very powerful UniProt-provided mapping file: its Ensembl mappings are a semicolon-delimeted list of Emsembl IDs containing
    # attached isoforms (and not all UniProtKB-AC identifiers have those!) so we'll need to do some extra post-processing.
    idmapping_selected_df = pandas.read_csv(
        interactome_folder / "HUMAN_9606_idmapping_selected.tsv",
        header=None,
        # See directory.py for the README associated with this mapping file.
        usecols=[0, 18],
        names=["UniProtKB-AC", "Ensembl"],
        sep="\t",
    )
    idmapping_selected_df = idmapping_selected_df[idmapping_selected_df["Ensembl"].notnull()]
    # Handle our ;-delimited list
    idmapping_selected_df['Ensembl'] = idmapping_selected_df['Ensembl'].str.split("; ")
    idmapping_selected_df = idmapping_selected_df.explode('Ensembl')
    # Drop isoforms
    idmapping_selected_df['Ensembl'] = idmapping_selected_df['Ensembl'].str.split('.').str[0]

    tf_df = tf_df.merge(idmapping_selected_df, on='Ensembl', how='inner')
    tf_df = tf_df.explode('UniProtKB-AC')
    tf_df = tf_df.fillna('NA')
    tf_df.to_csv(interactome_folder / "Homo_sapiens_TF_Uniprot.txt", header = True, sep = "\t", index = False)

if __name__ == "__main__":
    main()
