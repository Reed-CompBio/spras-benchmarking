from pathlib import Path
import pandas

egfr_directory = Path(__file__).parent.resolve() / ".."


def main():
    interactome_df = pandas.read_csv(egfr_directory / "raw" / "9606.protein.links.full.txt", sep=" ")
    interactome_df["protein1"] = interactome_df["protein1"].astype(str).str.removeprefix("9606.")
    interactome_df["protein2"] = interactome_df["protein2"].astype(str).str.removeprefix("9606.")
    # STRING provides two links (interactome) files: links and links.full. For convenience in other datasets,
    # we only download and use links.full, even though links would suffice. Due to that, we restrict to the
    # subset of headers in links before saving the interactome.
    interactome_df = interactome_df[["protein1", "protein2", "combined_score"]]
    interactome_df["Direction"] = "U"

    (egfr_directory / "processed").mkdir(exist_ok=True)
    interactome_df.to_csv(egfr_directory / "processed" / "ensp" / "interactome.tsv", index=False, header=False, sep="\t")


if __name__ == "__main__":
    main()
