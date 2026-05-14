from pathlib import Path
import pandas

from tools.normalize.interactome import deduplicate_edges

depmap_directory = Path(__file__).parent.resolve() / ".."


def main():
    (depmap_directory / "processed").mkdir(exist_ok=True)

    interactome_df = pandas.read_csv(depmap_directory / "raw" / "9606.protein.physical.links.full.txt", sep=" ")
    print(interactome_df)

    # rename the columns
    interactome_df = interactome_df.rename(columns={"protein1": "Interactor1", "protein2": "Interactor2", "combined_score": "Weight"})
    interactome_df["Interactor1"] = interactome_df["Interactor1"].astype(str).str.removeprefix("9606.")
    interactome_df["Interactor2"] = interactome_df["Interactor2"].astype(str).str.removeprefix("9606.")
    interactome_df = interactome_df[["Interactor1", "Interactor2", "Weight"]]
    interactome_df["Direction"] = "U"

    # deduplicate edges in the interactome
    interactome_df, _ = deduplicate_edges(interactome_df)
    interactome_df.sort_values('Weight',inplace=True)

    # save SPRAS formatted interactome
    interactome_df.to_csv(depmap_directory / "processed" / "interactome.tsv", index=False, header=False, sep="\t")


if __name__ == "__main__":
    main()


