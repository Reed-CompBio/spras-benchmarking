import pandas
from pathlib import Path

hiv_path = Path(__file__).parent.resolve().parent


def process_prizes(prizes: pandas.DataFrame):
    # Some proteins in the original prize files have the syntax `majorIdentifier-N` where N denotes isoforms.
    # We don't particurarly care about any particular isoform when doing pathway reconstruction,
    # so we treat -N isoforms as duplicates and remove them.
    prizes["Uniprot"] = prizes["Uniprot"].str.split("-", expand=False).str[0]

    # We sort for the highest Prize for all of the isoform (and non-isoform) variants
    # to make the output more readable.
    prizes = prizes.sort_values("Prize", ascending=False).drop_duplicates("Uniprot").sort_index()

    return prizes


def main():
    # Follow `Snakefile` or the README for information about these two files.
    prize_05 = process_prizes(pandas.read_csv(hiv_path / "raw" / "prize_05.tsv", sep="\t"))
    prize_060 = process_prizes(pandas.read_csv(hiv_path / "raw" / "prize_060.tsv", sep="\t"))

    prize_060_nodes = prize_060["Uniprot"].tolist()
    prize_05_nodes = prize_05["Uniprot"].tolist()
    node_set = list(set(prize_05_nodes + prize_060_nodes))

    # Save files to the intermediate path
    intermediate_path = hiv_path / "intermediate"
    intermediate_path.mkdir(exist_ok=True)
    prize_05.to_csv(intermediate_path / "prize_05.tsv", index=False, sep="\t")
    prize_060.to_csv(intermediate_path / "prize_060.tsv", index=False, sep="\t")
    (intermediate_path / "node_set.txt").write_text("\n".join(node_set))


if __name__ == "__main__":
    main()
