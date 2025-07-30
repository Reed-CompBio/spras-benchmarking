import pickle
from pathlib import Path
import os

current_directory = Path(os.path.dirname(os.path.realpath(__file__)))
processed_directory = current_directory.parent / "processed"

def main():
    with open("Pickles/UniprotIDs.pkl", "rb") as file:
        UniprotIDs = pickle.load(file)

    UIDs = UniprotIDs["UniprotIDs"]
    UMap = UniprotIDs["UniprotMap"]

    with open("Pickles/NodeIDs.pkl", "rb") as file2:
        prizes = pickle.load(file2)

    prize_05 = prizes["prize_05"]
    prize_060 = prizes["prize_060"]

    prize_05["Uniprot"] = prize_05["Uniprot"].apply(lambda x: UMap.get(x))
    prize_060["Uniprot"] = prize_060["Uniprot"].apply(lambda x: UMap.get(x))

    prize_05.columns = ["NODEID", "prize"]
    prize_060.columns = ["NODEID", "prize"]

    prize_05.to_csv(processed_directory / "processed_prize_05.txt", sep="\t", header=True, index=False)
    prize_060.to_csv(processed_directory / "processed_prize_060.txt", sep="\t", header=True, index=False)

if __name__ == '__main__':
    main()
