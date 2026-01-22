import pickle
from pathlib import Path

current_directory = Path(__file__).parent.resolve()
processed_directory = current_directory.parent / "processed"

def main():
    # See name_mapping.py for the origin of this file
    with open(current_directory / ".." / "pickles" / "UniprotIDs.pkl", "rb") as file:
        UniprotIDs = pickle.load(file)
    UMap = UniprotIDs["UniprotMap"]

    with open(current_directory / ".." / "pickles" / "NodeIDs.pkl", "rb") as file2:
        prizes = pickle.load(file2)
    prize_05 = prizes["prize_05"]
    prize_060 = prizes["prize_060"]

    prize_05["Uniprot"] = prize_05["Uniprot"].apply(lambda x: UMap.get(x))
    prize_060["Uniprot"] = prize_060["Uniprot"].apply(lambda x: UMap.get(x))

    # Format with SPRAS column names
    prize_05.columns = ["NODEID", "prize"]
    prize_060.columns = ["NODEID", "prize"]

    prize_05.to_csv(processed_directory / "processed_prizes_05.txt", sep="\t", header=True, index=False)
    prize_060.to_csv(processed_directory / "processed_prizes_060.txt", sep="\t", header=True, index=False)

if __name__ == '__main__':
    main()
