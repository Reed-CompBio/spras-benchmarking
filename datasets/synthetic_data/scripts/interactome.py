import pandas
from pathlib import Path
from tools.mapping.ensembl_uniprot import idmapping_as_ensp_uniprot_mapping, idmapping_uniprot_mapping

current_directory = Path(__file__).parent.resolve()
interactome_folder = current_directory / ".." / "raw" / "human-interactome"

def main():
    # Convert the interactome to SPRAS format
    print("Reading interactome...")
    interactome_df = pandas.read_csv(
        current_directory / ".." / "raw" / "9606.protein.links.full.v12.0.txt", sep=" ", usecols=["protein1", "protein2", "combined_score"]
    )
    interactome_df.columns = ["Protein1", "Protein2", "Weight"]

    # We also want to representatively remove a certain percentage of elements from the interactome,
    # to make sure our interactome downsampling preserves edge weight distributions
    # (we don't care to preserve other major topological properties just yet.)
    # since this file is large, we opt for streaming the interactome for removing edges instead

    print("Initially processing interactome...")
    interactome_df["Weight"] = interactome_df["Weight"].div(1000)  # scores are from 1-1000: we normalize from 0-1.
    interactome_df["Direction"] = "U"
    print("Sorting interactome...")
    interactome_df = interactome_df.sort_values("Weight", kind="stable")
    interactome_df = interactome_df.reset_index(drop=True)

    print("Fetching mapping data...")

    # Mapping ENSP IDs to ENSG IDs through the STRING aliases file
    string_aliases = pandas.read_csv(current_directory / ".." / "raw" / "9606.protein.aliases.txt", sep="\t", usecols=["#string_protein_id", "alias"])
    string_aliases.columns = ["ENSG", "ENSP"]
    string_aliases = string_aliases.drop_duplicates()

    # (ENSG) idmapping -> (ENSG <-> ENSP) -> (ENSP) idmapping
    idmapping_df = idmapping_uniprot_mapping(interactome_folder / "HUMAN_9606_idmapping_selected.tsv")
    idmapping_df = idmapping_as_ensp_uniprot_mapping(idmapping_df)

    print("Mapping interactome...")
    # We also use astype(str) as these are read as numpy objects for convenience, but this messes with merging
    interactome_df["Protein1"] = interactome_df["Protein1"].str.removeprefix("9606.").astype(str)
    interactome_df["Protein2"] = interactome_df["Protein2"].str.removeprefix("9606.").astype(str)

    interactome_df = interactome_df.merge(idmapping_df, left_on="Protein1", right_on="Ensembl", how="left") \
        .drop(columns=["Protein1", "Ensembl"]) \
        .rename(columns={"UniProtKB-AC": "Protein1"}) \
        .drop_duplicates()
    interactome_df = interactome_df.merge(idmapping_df, left_on="Protein2", right_on="Ensembl", how="left") \
        .drop(columns=["Protein2", "Ensembl"]) \
        .rename(columns={"UniProtKB-AC": "Protein2"})

    interactome_df = interactome_df.dropna(subset=["Protein1", "Protein2"]).reset_index(drop=True)
    interactome_df = interactome_df[["Protein1", "Protein2", "Weight", "Direction"]]

    print("Counting weight counts...")
    interactome_df["Weight"].value_counts(sort=False).to_csv(current_directory / ".." / "processed" / "weight-counts.tsv", sep="\t")

    print("Saving interactome...")
    interactome_df.to_csv(current_directory / ".." / "processed" / "interactome.tsv", sep="\t", header=False, index=False)

if __name__ == "__main__":
    main()
