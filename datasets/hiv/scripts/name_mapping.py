"""
Maps UniProtKB AC/ID to UniProtKB/Swiss-Prot.
"""

from pathlib import Path
import pandas

hiv_directory = Path(__file__).parent.resolve().parent

def main():
    # This is 1 of two major exceptions to this being example code from UniProt.
    # See prepare.py for the NodeIDs generation: this is the deduplicated list of node IDs
    # from the two prize files in `raw`.
    node_ids = (hiv_directory / "intermediate" / "node_set.txt").read_text().split("\n")

    idmapping = pandas.read_csv(hiv_directory / "raw" / "HUMAN_9606_idmapping.tsv",
                                sep='\t', header=None, names=["UniProtKB", "Type", "UniProtKB-ID"])
    idmapping = idmapping[idmapping["Type"] == "UniProtKB-ID"]
    idmapping = idmapping.drop(columns="Type")
    idmapping = idmapping[idmapping["UniProtKB"].isin(node_ids)]
    idmapping.to_csv(hiv_directory / "intermediate" / "mapping.tsv", index=False, sep='\t')

if __name__ == "__main__":
    main()
