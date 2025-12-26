import pandas as pd
import os
from pathlib import Path

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

    # Get the BioMart ENSP -> ENSG mapping
    biomart_data = pd.read_csv(diseases_path / "raw" / "ensg-ensp.tsv", sep="\t", names=["ENSP", "ENSG"])

    # The DISEASES data is in the ENSP namespace, but we want to work in ENSG.
    knowledge_mapped = text_mining.merge(biomart_data, left_on="geneID", right_on="ENSP", how="inner")
    text_mining_mapped = knowledge.merge(biomart_data, left_on="geneID", right_on="ENSP", how="inner")

    knowledge_mapped = (
        knowledge_mapped.sort_values("confidenceScore", ascending=False).drop_duplicates(subset=["ENSG", "diseaseID"], keep=False).sort_index()
    )
    text_mining_mapped = (
        text_mining_mapped.sort_values("confidenceScore", ascending=False).drop_duplicates(subset=["ENSG", "diseaseID"], keep=False).sort_index()
    )

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

    # Mapping ENSG IDs to STRING IDs through the STRING aliases file
    # given our ENSG and ENSP (non one-to-one!) mapping `string_aliases`,

    # NOTE: the STRING API call to map genes to proteins
    # also does text search, which brings up more false positives than true positives: because
    # of this, we specifically only care about ENSG -> ENSP and nothing greater.
    string_aliases = pd.read_csv(
        diseases_path / ".." / ".." / "databases" / "string" / "9606.protein.aliases.v12.0.txt", sep="\t", usecols=["#string_protein_id", "alias"]
    )
    string_aliases.columns = ["str_id", "ENSP"]
    string_aliases = string_aliases.drop_duplicates()

    GS_string_df = GS_combined_threshold.merge(string_aliases, on="ENSP", how="inner")
    GS_string_df = GS_string_df.drop_duplicates(subset=["ENSG", "ENSP", "geneName", "diseaseID", "diseaseName"])

    GS_string_df.to_csv(diseases_path / "data" / "gold_standard.csv", index=False)


if __name__ == "__main__":
    main()
