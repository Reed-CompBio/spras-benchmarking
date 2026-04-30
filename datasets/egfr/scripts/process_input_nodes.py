import pandas
from pathlib import Path

egfr_directory = Path(__file__).parent.resolve() / ".."

SOURCE = "EGF_HUMAN"


def main():
    prizes = pandas.read_csv(egfr_directory / "raw" / "egfr-prizes.txt", sep="\t", header=None, names=["NODEID", "prize"])
    prizes = prizes.loc[~prizes["NODEID"].str.endswith("_PSEUDONODE")]
    prizes["target"] = "True"
    # `SOURCE` is designated as the source node and assigned the highest score to anchor reconstruction at the perturbation.
    prizes = pandas.concat(
        [prizes, pandas.DataFrame({"NODEID": [SOURCE], "prize": [prizes.loc[:, "prize"].max()], "dummy": ["True"], "source": ["True"]})],
        ignore_index=True,
    )
    prizes["active"] = "True"

    prizes.to_csv(egfr_directory / "preprocessed" / "uniprot" / "input-nodes.txt", index=False, sep="\t")


if __name__ == "__main__":
    main()
