# This prepares prizes with dummy nodes.
import pandas as pd
from pathlib import Path

data_directory = Path(__file__).parent.resolve()

if __name__ == "__main__":
    # Get the raw prizes DF
    prizes = data_directory / "raw" / "prizes.txt"
    prizes_df = pd.read_csv(prizes, sep="\t", header=None, names=["NODEID", "prize"])

    # Use the manually curated prize info
    # TODO: where did this come from? Nothing in the original documents even mention this.
    prizes_df2 = pd.DataFrame(data={"NODEID": ["YGR014W", "YDR420W", "YER118C"], "prize": 10.051863}, index=[1596, 1597, 1598])

    new_prizes_path = data_directory / "processed" / "prizes1_dummies.txt"
    new_prizes = pd.concat([prizes_df, prizes_df2])
    new_prizes.to_csv(new_prizes_path, sep="\t", index=False, columns=["NODEID", "prize"], header=["NODEID", "prize"])
