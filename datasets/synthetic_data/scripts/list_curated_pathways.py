import json
from pathlib import Path
from jsonc_parser.parser import JsoncParser

from datasets.synthetic_data.util.parse_pc_pathways import parse_pc_pathways

synthetic_directory = Path(__file__).parent.parent.resolve()


def main():
    # TODO: pass as arguments
    pathways_df = parse_pc_pathways(synthetic_directory / "raw" / "pathways.txt")

    # We use the top-level pathways.jsonc, which is a hand-curated list of pathways, as it is not deterministically
    # automatable to decide whether or not a pathway is a signaling pathway. Yet.
    pathway_mapping: dict[str, str] = {}
    curated_pathways = JsoncParser.parse_file(synthetic_directory / "pathways.jsonc")
    for pathway in curated_pathways:
        selected_pathways = pathways_df.loc[pathways_df["DISPLAY_NAME"] == pathway].reset_index(drop=True)
        selected_pathways_count = len(selected_pathways.index)
        if selected_pathways_count != 1:
            raise RuntimeError(f"{pathway} references {selected_pathways_count} pathways, when we need to uniquely get one!")
        pathway_mapping[pathway] = selected_pathways["PATHWAY_URI"].loc[0]
    (synthetic_directory / "intermediate" / "curated_pathways_id_mapping.json").write_text(json.dumps(pathway_mapping, indent=4))


if __name__ == "__main__":
    main()
