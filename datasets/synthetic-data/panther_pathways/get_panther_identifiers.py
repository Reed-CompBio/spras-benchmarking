from io import StringIO
from pathlib import Path
import pandas

from jsonc_parser.parser import JsoncParser

current_directory = Path(__file__).parent.resolve()

def main():
    # We have two tables: the latter actually has more data, so we use that one instead.
    # These two tables are separated by two newlines.
    needle = "\n\n"
    _, text = (current_directory / 'raw' / 'pathways.txt').read_text().split(needle)

    pathways_df = pandas.read_csv(StringIO(text), sep='\t')
    pathways_df = pathways_df.loc[pathways_df["DATASOURCE"] == "PANTHER"]
    pathways_df = pathways_df.loc[pathways_df["NUM_DIRECT_COMPONENT_OR_STEP_PROCESSES"] != 0]
    pathways_df = pathways_df.reset_index(drop=True)
    pathways_df = pathways_df[["PATHWAY_URI", "DISPLAY_NAME"]]
    (current_directory / "intermediate").mkdir(exist_ok=True)
    pathways_df.to_csv(current_directory / "intermediate" / "all_pathways.tsv", index=False, sep='\t')

    # We use the top-level pathways.jsonc, which is a hand-curated list of pathways, as it is not deterministically
    # automatable to decide whether or not a pathway is a signaling pathway. Yet.
    pathway_mapping: dict[str, str] = {}
    curated_pathways = JsoncParser.parse_file(current_directory / "pathways.jsonc")
    for pathway in curated_pathways:
        selected_pathways = pathways_df.loc[pathways_df["DISPLAY_NAME"] == pathway].reset_index(drop=True)
        selected_pathways_count = len(selected_pathways.index)
        if selected_pathways_count != 1:
            raise RuntimeError(f"{pathway} references {selected_pathways_count} pathways, when we need to uniquely get one!")
        pathway_mapping[pathway] = selected_pathways["PATHWAY_URI"].loc[0]
    curated_pathway_df = pandas.DataFrame(pathway_mapping.items())
    curated_pathway_df.columns = ["Name", "ID"]
    curated_pathway_df.to_csv(current_directory / "intermediate" / "curated_pathways.tsv", index=False, sep='\t')

if __name__ == "__main__":
    main()
