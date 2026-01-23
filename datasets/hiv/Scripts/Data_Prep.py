import pandas
import pickle
import os

prize_05 = pandas.read_csv("raw/prize_05.csv", sep="\t", lineterminator="\n")
prize_060 = pandas.read_csv("raw/prize_060.csv", sep="\t", lineterminator="\n")

prize_05["Uniprot"] = prize_05["Uniprot"].str.split("-", expand=False).str[0]
prize_05 = prize_05.sort_values("Prize", ascending=False).drop_duplicates("Uniprot").sort_index()

prize_060["Uniprot"] = prize_060["Uniprot"].str.split("-", expand=False).str[0]
prize_060 = prize_060.sort_values("Prize", ascending=False).drop_duplicates("Uniprot").sort_index()

prize_060_nodes = prize_060["Uniprot"].tolist()
prize_05_nodes = prize_05["Uniprot"].tolist()

nodeset = list(set(prize_05_nodes + prize_060_nodes))

df = {"NodeIDs": nodeset, "prize_05": prize_05, "prize_060": prize_060}

if not os.path.exists("./Pickles"):
    os.makedirs("./Pickles")

with open("Pickles/NodeIDs.pkl", "wb") as file:
    pickle.dump(df, file)
