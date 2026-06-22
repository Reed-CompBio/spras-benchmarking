from pathlib import Path
import pandas

from tools.normalize.interactome import deduplicate_edges

egfr_directory = Path(__file__).parent.resolve() / ".."


def main():

    interactome_df = pandas.read_csv(egfr_directory / "raw" / "9606.protein.physical.links.txt", sep=" ")
    # Rename the columns both to stylistically keep it in-line with SPRAS and functionally for `deduplicate_edges`.
    interactome_df = interactome_df.rename(columns={"protein1": "Interactor1", "protein2": "Interactor2", "combined_score": "Weight"})
    # keeping IDs as ESNP
    interactome_df["Interactor1"] = interactome_df["Interactor1"].astype(str).str.removeprefix("9606.")
    interactome_df["Interactor2"] = interactome_df["Interactor2"].astype(str).str.removeprefix("9606.")
    interactome_df = interactome_df[["Interactor1", "Interactor2", "Weight"]]
    interactome_df["Direction"] = "U"
    # TODO: update the weights to be between 0 - 1 instead of 150 - 999 by / 1000? Or should we have this be specific to algorithms fixing this? 
    # For now I am converting it to be 0 - 1
    interactome_df["Weight"] = interactome_df["Weight"] / 1000


    # We normalize the interactome (any final post-processing steps wanted/needed by SPRAS).
    (egfr_directory / "preprocessed").mkdir(exist_ok=True)
    interactome_df, _ = deduplicate_edges(interactome_df)
    interactome_df.to_csv(egfr_directory / "preprocessed" / "ensp" / "interactome.tsv", index=False, header=False, sep="\t")


if __name__ == "__main__":
    main()
