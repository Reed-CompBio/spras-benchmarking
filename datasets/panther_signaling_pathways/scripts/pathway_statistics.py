import sys
from pathlib import Path
import pandas

panther_directory = Path(__file__).parent.parent.resolve()
intermediate_directory = panther_directory / "intermediate"
raw_directory = panther_directory / "raw"
processed_directory = panther_directory / "processed"


def process_pathway(pathway: str) -> dict:
    pathway_folder = intermediate_directory / pathway
    pathway_file = pandas.read_csv(pathway_folder / "pathway.csv", sep="\t")
    input_nodes = pandas.read_csv(pathway_folder / "input_nodes.csv", sep="\t")

    sources = list(input_nodes[input_nodes["sources"] == True]["NODEID"])
    targets = list(input_nodes[input_nodes["targets"] == True]["NODEID"])

    n_sources = len(sources)
    n_targets = len(targets)

    return {
        "Name": pathway,
        "Len_pathway": len(pathway_file),
        "Sources": n_sources,
        "Targets": n_targets,
        "Ratio Sources/Targets": (float(n_sources) / float(n_targets)) if n_targets != 0 else float("nan"),
        "Ratio Targets/Sources": (float(n_targets) / float(n_sources)) if n_sources != 0 else float("nan"),
    }


def main():
    pathways = sys.argv[1:]
    if not pathways:
        print("Usage: pathway_statistics.py <pathway1> <pathway2> ...")
        sys.exit(1)

    data_entries = [process_pathway(p) for p in pathways]

    data_df = pandas.DataFrame(data_entries)
    data_df.to_csv(processed_directory / "full_stats.tsv", sep="\t", index=False)

    filtered_df = (
        data_df
        .loc[data_df["Sources"] != 0]
        .loc[data_df["Targets"] != 0]
        .reset_index(drop=True)
    )
    filtered_df.to_csv(processed_directory / "filtered_stats.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
