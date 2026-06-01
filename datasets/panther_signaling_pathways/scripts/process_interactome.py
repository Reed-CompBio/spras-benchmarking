# get each of the interactomes
import pandas as pd
from pathlib import Path

from tools.normalize.interactome import deduplicate_edges

panther_directory = Path(__file__).parent.parent.resolve()
processed_directory = panther_directory / "processed"
raw_directory = panther_directory / "raw"

def main():
    (processed_directory).mkdir(exist_ok=True)

    interactome_df = pd.read_csv(
        raw_directory / "9606.protein.links.full.v12.0.txt", sep=" ", usecols=["protein1", "protein2", "combined_score"]
    )

    interactome_df = interactome_df.rename(columns={"protein1": "Interactor1", "protein2": "Interactor2", "combined_score": "Weight"})
    interactome_df = interactome_df[interactome_df["Weight"] > 1]
    interactome_df["Weight"] = interactome_df["Weight"].div(1000) # combined_scores are between 1-1000: we normalize between 0-1
    interactome_df["Interactor1"] = interactome_df["Interactor1"].astype(str).str.removeprefix("9606.")
    interactome_df["Interactor2"] = interactome_df["Interactor2"].astype(str).str.removeprefix("9606.")
    interactome_df["Direction"] = "U"

    # deduplicate edges in the interactome
    interactome_df, _ = deduplicate_edges(interactome_df)
    interactome_df.sort_values('Weight',inplace=True)

    # save SPRAS formatted interactome
    interactome_df.to_csv(processed_directory / "interactome.tsv", index=False, header=False, sep="\t")

    # save weight count distributions
    interactome_df["Weight"].value_counts(sort=False).to_csv(processed_directory / "weight-counts.tsv", sep="\t")

if __name__ == "__main__":
    main()