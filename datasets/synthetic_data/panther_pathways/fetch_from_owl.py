from pathlib import Path
from paxtools.fetch import fetch
from datasets.synthetic_data.util.parse_pc_pathways import parse_pc_pathways

current_directory = Path(__file__).parent.resolve()


def main():
    pathways_df = parse_pc_pathways(current_directory / "raw" / "pathways.txt")
    print("Fetching pathways... [This may take some time. On the author's desktop machine, it took 15 minutes.]")
    (current_directory / "output").mkdir(exist_ok=True)
    fetch(
        current_directory / "raw" / "pc-biopax.owl",
        output=(current_directory / "output" / "pc-panther-biopax.owl"),
        uris=list(pathways_df["PATHWAY_URI"]),
        memory=f"{2 ** (16 - 1)}m",  # this is why we don't run this in CI! This is 32gb of memory.
    )


if __name__ == "__main__":
    main()
