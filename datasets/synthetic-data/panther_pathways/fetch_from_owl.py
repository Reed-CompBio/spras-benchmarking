from pathlib import Path
from paxtools.fetch import fetch
from ..util.parse_pc_pathways import parse_pc_pathways

current_directory = Path(__file__).parent.resolve()

def main():
    pathways_df = parse_pc_pathways(current_directory / 'raw' / 'pathways.txt')
    print(list(pathways_df[["PATHWAY_URI"]]))
    # fetch(current_directory / 'raw' / 'pc-biopax.owl', output="pc-biopax-selected.owl")
    pass

if __name__ == "__main__":
    main()