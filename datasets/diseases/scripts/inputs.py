import pandas as pd
import requests
import pickle


def main():
    tiga = pd.read_csv("datasets/diseases/raw/tiga_gene-trait_stats.tsv", sep="\t")
    tiga = tiga[["ensemblId", "efoId", "trait", "n_snp", "n_snpw"]]
    tiga = tiga.drop_duplicates(subset=["ensemblId", "trait"])

    human_do = pd.read_csv("datasets/diseases/raw/HumanDO.tsv", sep="\t")
    human_do = human_do.drop_duplicates(subset="label")
    human_do = human_do[["id", "label"]]

    tiga_do = tiga.merge(human_do, left_on="trait", right_on="label", how="inner", validate="m:1")

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

    with open("datasets/diseases/pickles/inputs.pkl", "wb") as file:
        pickle.dump(tiga_string_df, file)

if __name__ == "__main__":
    main()
