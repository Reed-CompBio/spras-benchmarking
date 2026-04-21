from pathlib import Path
import pandas

from tools.interactome import normalize_interactome

egfr_directory = Path(__file__).parent.resolve() / ".."


def main():
    interactome_df = pandas.read_csv(egfr_directory / "raw" / "9606.protein.links.full.txt", sep=" ")
    # Rename the columns both to stylistically keep it in-line with SPRAS and functionally for `normalize_interactome`.
    interactome_df = interactome_df.rename(columns={"protein1": "Interactor1", "protein2": "Interactor2", "combined_score": "Weight"})
    interactome_df["Interactor1"] = interactome_df["Interactor1"].astype(str).str.removeprefix("9606.")
    interactome_df["Interactor2"] = interactome_df["Interactor2"].astype(str).str.removeprefix("9606.")
    # STRING provides two links (interactome) files: links and links.full. For convenience in other datasets,
    # we only download and use links.full, even though links would suffice. Due to that, we restrict to the
    # subset of headers in links before saving the interactome.
    interactome_df = interactome_df[["Interactor1", "Interactor2", "Weight"]]
    interactome_df["Direction"] = "U"

    # We normalize the interactome (any final post-processing steps wanted/needed by SPRAS).
    (egfr_directory / "processed").mkdir(exist_ok=True)
    interactome_df, _ = normalize_interactome(interactome_df)
    interactome_df.to_csv(egfr_directory / "processed" / "ensp" / "interactome.tsv", index=False, header=False, sep="\t")


if __name__ == "__main__":
    main()
