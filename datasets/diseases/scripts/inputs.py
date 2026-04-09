from pathlib import Path
import pandas as pd

dir_path = Path(__file__).parent.resolve()

diseases_path = Path(dir_path, "..")
(diseases_path / "data").mkdir(exist_ok=True, parents=True)

'''
Generates input node set files based on TIGA trait-gene associations and Disease Ontology
annotations. 
TODO: not all of these input node set files are necessary for benchmarking; only those with
disease validation sets. That is processed in `files.py`. 
'''
def main():
    # See ../Snakefile for information on this file's origin
    tiga = pd.read_csv(diseases_path / "raw" / "tiga_gene-trait_stats.tsv", sep="\t")
    tiga = tiga[["ensemblId", "trait", "n_snp", "n_snpw"]]
    tiga = tiga.drop_duplicates(subset=["ensemblId", "trait"])

    # See ../Snakefile for information on this file's origin
    human_do = pd.read_csv(diseases_path / "raw" / "HumanDO.tsv", sep="\t")
    human_do = human_do.drop_duplicates(subset="label")
    human_do = human_do[["id", "label"]]

    # Get all TIGA trait-gene associations where the traits correspond to Disease Ontology (DO) terms.
    tiga_do = tiga.merge(human_do, left_on="trait", right_on="label", how="inner", validate="many_to_one")

    # Mapping ENSG IDs to ENSP IDs through the STRING aliases file
    # given our ENSG and ENSP (non one-to-one!) mapping `string_aliases`.
    string_aliases = pd.read_csv(diseases_path / "raw" / "9606.protein.aliases.txt", sep="\t", usecols=["#string_protein_id", "alias"])
    string_aliases.columns = ["str_id", "ENSP"]
    string_aliases = string_aliases.drop_duplicates()

    # We can create our TIGA-mapped file through merging on the ENSPs
    tiga_string_df = tiga_do.merge(string_aliases, left_on="ensemblId", right_on="ENSP", how="inner")

    ## Write file to inputs.csv
    tiga_string_df.to_csv(diseases_path / "data" / "inputs.csv", index=False)


if __name__ == "__main__":
    main()
