import pandas
from pathlib import Path
from ensembl_uniprot_mapping import idmapping_uniprot_mapping, idmapping_as_ensg_uniprot_mapping

current_directory = Path(__file__).parent.resolve()

interactome_folder = current_directory / ".." / "raw" / "human-interactome"


def main():
    tf_df = pandas.read_csv(interactome_folder / "Homo_sapiens_TF.tsv", sep="\t", header=0)
    idmapping_selected_df = idmapping_uniprot_mapping(interactome_folder / "HUMAN_9606_idmapping_selected.tsv")
    idmapping_selected_df = idmapping_as_ensg_uniprot_mapping(idmapping_selected_df)
    tf_df = tf_df.merge(idmapping_selected_df, on="Ensembl", how="inner")
    tf_df = tf_df.explode("UniProtKB-AC")
    tf_df = tf_df.fillna("NA")
    tf_df.to_csv(interactome_folder / "Homo_sapiens_TF_Uniprot.tsv", header=True, sep="\t", index=False)


if __name__ == "__main__":
    main()
