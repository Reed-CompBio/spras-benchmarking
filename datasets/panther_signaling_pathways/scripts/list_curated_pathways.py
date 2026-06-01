import json
from pathlib import Path
from os import PathLike
from jsonc_parser.parser import JsoncParser
from io import StringIO
import pandas

panther_directory = Path(__file__).parent.parent.resolve()

def parse_pc_pathways(pathways_txt_path: str | PathLike) -> pandas.DataFrame:
    """
    Parses the pathways.txt file from the PathwayCommons FTP store into a Dataframe of the sort
    PATHWAY_URI\tDISPLAY_NAME

    such that DATASOURCE is filtered by PANTHER.
    """
    # We have two tables: the latter actually has more data, so we use that one instead.
    # These two tables are separated by two newlines.
    needle = "\n\n"
    _, text = Path(pathways_txt_path).read_text().split(needle)

    pathways_df = pandas.read_csv(StringIO(text), sep="\t")
    pathways_df = pathways_df.loc[pathways_df["DATASOURCE"] == "PANTHER"]
    pathways_df = pathways_df.loc[pathways_df["NUM_DIRECT_COMPONENT_OR_STEP_PROCESSES"] != 0]
    pathways_df = pathways_df.reset_index(drop=True)
    pathways_df = pathways_df[["PATHWAY_URI", "DISPLAY_NAME"]]
    return pathways_df

def main():

    # this file is pathway hierarchy table from Pathway Commons of metadata is shipped with release v14
    # aggregates the pathways from their sources and assigns each a URL
    pathways_df = parse_pc_pathways(panther_directory / "raw" / "pathways.txt")

    # We use the top-level pathways.jsonc, which is a hand-curated list of signaling pathways,
    # Not automatable to decide whether or not a pathway is a signaling pathway from Panther Pathways database
    pathway_mapping: dict[str, str] = {}
    curated_pathways = JsoncParser.parse_file(panther_directory / "pathways.jsonc")
    for pathway in curated_pathways:
        selected_pathways = pathways_df.loc[pathways_df["DISPLAY_NAME"] == pathway].reset_index(drop=True)
        selected_pathways_count = len(selected_pathways.index)
        if selected_pathways_count != 1:
            raise RuntimeError(f"{pathway} references {selected_pathways_count} pathways, when we need to uniquely get one!")
        pathway_mapping[pathway] = selected_pathways["PATHWAY_URI"].loc[0]

    (panther_directory / "intermediate" / "curated_pathways_id_mapping.json").write_text(json.dumps(pathway_mapping, indent=4))

if __name__ == "__main__":
    main()
