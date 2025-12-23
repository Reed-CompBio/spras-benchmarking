from pathlib import Path
import pandas as pd
import requests
import os

# https://stackoverflow.com/a/5137509/7589775
dir_path = os.path.dirname(os.path.realpath(__file__))

diseases_path = Path(dir_path, "..")
(diseases_path / "data").mkdir(exist_ok=True, parents=True)


def main():
    # See fetch.py for information on this file's origin
    tiga = pd.read_csv(diseases_path / "raw" / "tiga_gene-trait_stats.tsv", sep="\t")
    tiga = tiga[["ensemblId", "trait", "n_snp", "n_snpw"]]
    tiga = tiga.drop_duplicates(subset=["ensemblId", "trait"])

    # See fetch.py for information on this file's origin
    human_do = pd.read_csv(diseases_path / "raw" / "HumanDO.tsv", sep="\t")
    human_do = human_do.drop_duplicates(subset="label")
    human_do = human_do[["id", "label"]]

    tiga_do = tiga.merge(human_do, left_on="trait", right_on="label", how="inner", validate="many_to_one")

    # Mapping ENSG IDs to STRING IDs through the STRING aliases file
    # given our ENSG and ENSP (non one-to-one!) mapping `string_aliases`,
    string_aliases = pd.read_csv(
        diseases_path / ".." / ".." / "databases" / "string" / "9606.protein.aliases.v12.0.txt",
        sep="\t", usecols=["#string_protein_id", "alias"])
    string_aliases.columns = ["ensp", "ensg"]
    string_aliases = string_aliases.drop_duplicates()
    string_aliases.to_csv(diseases_path / 'test.csv')

    # We can create our TIGA-mapped file through merging on the ENSPs
    tiga_string_df = tiga_do.merge(string_aliases, left_on="ensemblId", right_on="ensg", how='inner')
    tiga_string_df.to_csv(diseases_path / "data" / "inputs.csv", index=False)

if __name__ == "__main__":
    main()
