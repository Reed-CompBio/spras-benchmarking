import pandas as pd
import requests
import os
from pathlib import Path


def gprofiler_convert(ids: list, namespace: str, df: pd.DataFrame) -> pd.DataFrame:
    """
    Converts a list of IDs to a namespace, and associates it with a dataframe
    (see code for more details) using the gprofiler API.
    """

    # See the associated API documentation at https://biit.cs.ut.ee/gprofiler/page/apis.
    r = requests.post(
        url="https://biit.cs.ut.ee/gprofiler/api/convert/convert/",
        json={
            "organism": "hsapiens",
            "target": namespace,
            "query": ids,
        },
    )
    results = r.json()["result"]
    mapping = {}
    for x in results:
        if x["converted"] != "None":
            mapping.update({x["incoming"]: x["converted"]})

    mapping_df = pd.DataFrame.from_dict(mapping.items())
    mapping_df.columns = ["ENSP", "ENSG"]
    output_df = df.merge(mapping_df, left_on="geneID", right_on="ENSP", how="inner", validate="m:1")
    output_df = output_df.sort_values("confidenceScore", ascending=False).drop_duplicates(subset=["ENSG", "diseaseID"], keep=False).sort_index()

    return output_df


# https://stackoverflow.com/a/5137509/7589775
dir_path = os.path.dirname(os.path.realpath(__file__))

diseases_path = Path(dir_path, "..")


def main():
    # Get our data from `fetch.py`
    text_mining = pd.read_csv(diseases_path / "raw" / "human_disease_textmining_filtered.tsv", sep="\t")
    knowledge = pd.read_csv(diseases_path / "raw" / "human_disease_knowledge_filtered.tsv", sep="\t")
    # and correctly map their columns
    text_mining.columns = ["geneID", "geneName", "diseaseID", "diseaseName", "zScore", "confidenceScore", "sourceUrl"]
    knowledge.columns = ["geneID", "geneName", "diseaseID", "diseaseName", "sourceDB", "evidenceType", "confidenceScore"]

    # The DISEASES data is in the ENSP namespace, but we want to work in ENSG.
    knowledge_mapped = gprofiler_convert(list(text_mining["geneID"]), "ENSG", text_mining)
    text_mining_mapped = gprofiler_convert(list(knowledge["geneID"]), "ENSG", knowledge)

    inner = text_mining_mapped.merge(knowledge_mapped, on=["ENSG", "diseaseID"], how="inner")
    inner["confidenceScore"] = inner.apply(lambda x: max(x.confidenceScore_x, x.confidenceScore_y), axis=1)
    inner = inner.rename(columns={"ENSP_x": "ENSP", "geneName_x": "geneName", "diseaseName_x": "diseaseName", "geneID_x": "geneID"})
    inner = inner[["ENSG", "ENSP", "geneName", "diseaseID", "diseaseName", "confidenceScore"]]

    txt_only = text_mining_mapped.merge(knowledge_mapped, on=["ENSG", "diseaseID"], how="left")
    txt_only = txt_only[txt_only["confidenceScore_y"].isna()]
    txt_only = txt_only.rename(
        columns={"confidenceScore_x": "confidenceScore", "ENSP_x": "ENSP", "geneName_x": "geneName", "diseaseName_x": "diseaseName"}
    )
    txt_only = txt_only[["ENSG", "ENSP", "geneName", "diseaseID", "diseaseName", "confidenceScore"]]

    kn_only = text_mining_mapped.merge(knowledge_mapped, on=["ENSG", "diseaseID"], how="right")
    kn_only = kn_only[kn_only["confidenceScore_x"].isna()]
    kn_only = kn_only.rename(
        columns={"confidenceScore_y": "confidenceScore", "ENSP_y": "ENSP", "geneName_y": "geneName", "diseaseName_y": "diseaseName"}
    )
    kn_only = kn_only[["ENSG", "ENSP", "geneName", "diseaseID", "diseaseName", "confidenceScore"]]

    df_list = [inner, txt_only, kn_only]
    GS_ids = pd.concat(df_list)

    GS_score_threshold = GS_ids.loc[(GS_ids["confidenceScore"] >= 4)]
    GS_score_group = GS_ids.groupby("diseaseName")
    GS_score_dict = {k: v for k, v in GS_score_group}
    GS_score_count = {x: len(GS_score_dict[x]) for x in GS_score_dict.keys()}
    GS_count_threshold = {k: v for (k, v) in GS_score_count.items() if (v > 10)}
    GS_combined_threshold = GS_score_threshold.loc[GS_score_threshold["diseaseName"].isin(list(GS_count_threshold.keys()))]

    string_api_url = "https://version-12-0.string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"
    str_params = {
        "identifiers": "\r".join(list(GS_combined_threshold["ENSP"])),
        "species": 9606,
        "echo_query": 1,
    }
    request_url = "/".join([string_api_url, output_format, method])
    string_results = requests.post(request_url, data=str_params)
    string_results = string_results.text.strip()

    (diseases_path / "data").mkdir(exist_ok=True)

    # Archive the STRING data
    (diseases_path / "data" / "string-call.txt").write_text(string_results)
    string_map = {}
    for line in string_results.split("\n"):
        l = line.split("\t")
        string_map.update({l[0]: l[2]})
    string_df = pd.DataFrame.from_dict(string_map.items())
    string_df.columns = ["ENSP", "str_id"]
    GS_string_df = GS_combined_threshold.merge(string_df, on="ENSP", how="inner")
    GS_string_df.to_csv(diseases_path / "data" / "gold_standard.csv", index=False)


if __name__ == "__main__":
    main()
