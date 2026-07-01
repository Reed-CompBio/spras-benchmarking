from pathlib import Path
import pandas
from tools.normalize.interactome import deduplicate_edges

diseases_path = Path(__file__).parent.parent.resolve()


def main():
    # See /cache/directory.py for information on how this was grabbed.
    interactome_df = pandas.read_csv(diseases_path / "raw" / "9606.protein.links.full.txt", sep=" ")
    interactome_df = interactome_df.rename(columns={"protein1": "Interactor1", "protein2": "Interactor2", "combined_score": "Weight"})
    # keeping IDs as ESNP
    interactome_df["Interactor1"] = interactome_df["Interactor1"].astype(str).str.removeprefix("9606.")
    interactome_df["Interactor2"] = interactome_df["Interactor2"].astype(str).str.removeprefix("9606.")
    interactome_df = interactome_df[["Interactor1", "Interactor2", "Weight"]]
    interactome_df["Direction"] = "U"
    # TODO: update the weights to be between 0 - 1 instead of 1 - 999 by / 1000? Or should we have this be specific to algorithms fixing this?
    # For now I am converting it to be 0 - 1
    interactome_df["Weight"] = interactome_df["Weight"] / 1000
    interactome_df,_ = deduplicate_edges(interactome_df)

    (diseases_path / "processed").mkdir(exist_ok=True)
    interactome_df.to_csv(diseases_path / "processed" / "string_interactome.tsv", sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
