from pathlib import Path
from jsonc_parser.parser import JsoncParser
import pandas

current_directory = Path(__file__).parent.resolve()

def main():
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
