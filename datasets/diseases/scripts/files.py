import pickle
import pandas as pd
from pathlib import Path
import os

# https://stackoverflow.com/a/5137509/7589775
dir_path = os.path.dirname(os.path.realpath(__file__))

diseases_path = Path(dir_path, '..')
(diseases_path / 'prize_files').mkdir(exist_ok=True, parents=True)
(diseases_path / 'GS_files').mkdir(exist_ok=True, parents=True)

def main():
    with open(diseases_path / "pickles" / "gold_standard.pkl", "rb") as file:
        # See gold_standard.py
        GS_string_df = pickle.load(file)

    with open(diseases_path / "pickles" / "inputs.pkl", "rb") as file:
        # See inputs.py
        tiga_string_df = pickle.load(file)

    GS_string_df = GS_string_df[GS_string_df["diseaseID"].isin(tiga_string_df["id"])]
    GS_combined_group = GS_string_df.groupby("diseaseName")
    GS_combined_dict = {k: v for k, v in GS_combined_group}

    tiga_filtered = tiga_string_df[tiga_string_df["id"].isin(GS_string_df["diseaseID"])]
    tiga_group = tiga_filtered.groupby("trait")
    tiga_dict = {k: v for k, v in tiga_group}
    tiga_count = {x: len(tiga_dict[x]) for x in tiga_dict.keys()}
    tiga_count_threshold = {k: v for (k, v) in tiga_count.items() if (v > 10)}

    tiga_threshold = tiga_filtered.loc[tiga_filtered["trait"].isin(list(tiga_count_threshold.keys()))]

    tiga_prizes = tiga_threshold.groupby("trait")
    tiga_prize_dict = {k: v for k, v in tiga_prizes}

    for disease in tiga_prize_dict.keys():
        df = tiga_prize_dict[disease]
        df = df[["str_id", "n_snpw"]]
        df = df.rename(columns={"str_id": "NODEID", "n_snpw": "prize"})
        df.to_csv(diseases_path / "prize_files" / f"{disease.replace(' ', '_')}_prizes.txt", sep="\t", index=False)

    for disease in GS_combined_dict.keys():
        df = GS_combined_dict[disease]
        df = df[["str_id"]]
        df.to_csv(diseases_path / "GS_files" / f"{disease.replace(' ', '_')}_GS.txt", sep="\t", index=False, header=None)

    # See /databases/stringdb.py for information on how this was grabbed.
    # 9606 is the organism code for homo sapiens and the required background interactome of DISEASES.
    string = pd.read_csv(
        diseases_path / '..' / '..' / 'databases' / 'string' / '9606.protein.links.v12.0.txt',
        sep=" ", skiprows=[0], header=None)
    # Threshold anything above a confidence score of 900 to trim down the background interactome
    string = string[string.iloc[:, 2] > 900]
    string = string.iloc[:, [0, 1]]
    string[len(string.columns)] = 1
    string.to_csv(diseases_path / "raw" / "string_interactome.txt", sep="\t", index=False, header=None)

if __name__ == "__main__":
    main()
