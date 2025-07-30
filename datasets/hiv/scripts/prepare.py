import pandas
from pathlib import Path
import pickle
import os

# https://stackoverflow.com/a/5137509/7589775
hiv_path = Path(os.path.dirname(os.path.realpath(__file__)), '..')

def main():
    prize_05 = pandas.read_csv(hiv_path / "raw" / "prize_05.csv", sep="\t", lineterminator="\n")
    prize_060 = pandas.read_csv(hiv_path / "raw" / "prize_060.csv", sep="\t", lineterminator="\n")

    prize_05["Uniprot"] = prize_05["Uniprot"].str.split("-", expand=False).str[0]
    prize_05 = prize_05.sort_values("Prize", ascending=False).drop_duplicates("Uniprot").sort_index()

    prize_060["Uniprot"] = prize_060["Uniprot"].str.split("-", expand=False).str[0]
    prize_060 = prize_060.sort_values("Prize", ascending=False).drop_duplicates("Uniprot").sort_index()

    prize_060_nodes = prize_060["Uniprot"].tolist()
    prize_05_nodes = prize_05["Uniprot"].tolist()

    nodeset = list(set(prize_05_nodes + prize_060_nodes))

    df = {"NodeIDs": nodeset, "prize_05": prize_05, "prize_060": prize_060}

    (hiv_path / "Pickles").mkdir(exist_ok=True)

    with open(hiv_path / "Pickles" / "NodeIDs.pkl", "wb") as file:
        pickle.dump(df, file)

if __name__ == '__main__':
    main()

