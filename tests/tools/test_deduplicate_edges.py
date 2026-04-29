import pandas
import numpy
from tools.normalize.interactome import normalize_interactome


def test_normalize_interactome():
    df = pandas.DataFrame(
        {"Interactor1": ["a", "b", "e", "d"], "Interactor2": ["b", "a", "d", "e"], "Weight": [1, 2, 3, 4], "Direction": ["U", "U", "D", "D"]}
    )

    assert len(df.index) == 4
    normalized_df, _ = normalize_interactome(df)
    assert len(normalized_df.index) == 3
    assert list(normalized_df.iloc[0]) == ["a", "b", numpy.int64(1), "U"]
