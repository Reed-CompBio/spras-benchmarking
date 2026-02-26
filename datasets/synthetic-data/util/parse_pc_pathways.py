from io import StringIO
from os import PathLike
from pathlib import Path

import pandas


def parse_pc_pathways(pathways_txt_path: str | PathLike) -> pandas.DataFrame:
    """
    Parses a pathways.txt file from the PathwayCommons FTP store into a Dataframe of the sort
    PATHWAY_URI\tDISPLAY_NAME

    such that DATASOURCE is filtered by PANTHER. TODO: generalize to other datasources
    """
    # We have two tables: the latter actually has more data, so we use that one instead.
    # These two tables are separated by two newlines.
    needle = "\n\n"
    _, text = Path(pathways_txt_path).read_text().split(needle)

    pathways_df = pandas.read_csv(StringIO(text), sep='\t')
    pathways_df = pathways_df.loc[pathways_df["DATASOURCE"] == "PANTHER"]
    pathways_df = pathways_df.loc[pathways_df["NUM_DIRECT_COMPONENT_OR_STEP_PROCESSES"] != 0]
    pathways_df = pathways_df.reset_index(drop=True)
    pathways_df = pathways_df[["PATHWAY_URI", "DISPLAY_NAME"]]
    return pathways_df
