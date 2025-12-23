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

    # https://string-db.org/cgi/help.pl?subpage=api%23mapping-identifiers
    # We need to map ENSG IDs to STRING IDs
    print("Fetching STRING IDs...")
    string_api_url = "https://version-12-0.string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"
    str_params = {
        "identifiers": "\r".join(list(tiga_do["ensemblId"])),
        "species": 9606,
        "echo_query": 1,
    }
    request_url = "/".join([string_api_url, output_format, method])
    string_results = requests.post(request_url, data=str_params)

    string_map: dict[str, str] = {}
    for line in string_results.text.strip().split("\n"):
        l = line.split("\t")
        string_map.update({l[0]: l[2]})
    string_df = pd.DataFrame.from_dict(string_map.items())
    string_df.columns = ["ENSP", "str_id"]

    tiga_string_df = tiga_do.merge(string_df, left_on="ensemblId", right_on="ENSP", how="inner")

    tiga_string_df.to_csv(diseases_path / "data" / "inputs.csv", index=False)


if __name__ == "__main__":
    main()
