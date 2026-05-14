import pandas
import numpy
from tools.normalize.interactome import deduplicate_edges


def test_deduplicate_edges():
    df = pandas.DataFrame(
        {"Interactor1": ["a", "b", "e", "d"], "Interactor2": ["b", "a", "d", "e"], "Weight": [1, 2, 3, 4], "Direction": ["U", "U", "D", "D"]}
    )

    print(df)
    assert len(df.index) == 4
    deduplicated_df, _ = deduplicate_edges(df)

    print(deduplicated_df)
    assert len(deduplicated_df.index) == 3
    assert list(deduplicated_df.iloc[0]) == ["d", "e", numpy.int64(4), "D"]
