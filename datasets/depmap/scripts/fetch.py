"""
Fetches the latest DepMap data we need

Download page: https://depmap.org/portal/data_page/?tab=allData
"""

from pathlib import Path
import os
from cache.directory import get_cache_item
from databases.util import uncompress

# https://stackoverflow.com/a/5137509/7589775
dir_path = os.path.dirname(os.path.realpath(__file__))

raw_dir = Path(dir_path, "..", "raw")


def main():
    raw_dir.mkdir(exist_ok=True)

    print("Fetching DepMap omics metadata")
    get_cache_item(["DepMap", "OmicsProfiles.csv"]).download(raw_dir / "OmicsProfiles.csv")

    print("Fetching DepMap gene dependency probability estimates...")
    get_cache_item(["DepMap", "CRISPRGeneDependency.csv"]).download(raw_dir / "CRISPRGeneDependency.csv")

    print("Fetching DepMap genotyped matrix...")
    get_cache_item(["DepMap", "OmicsSomaticMutationsMatrixDamaging.csv"]).download(raw_dir / "OmicsSomaticMutationsMatrixDamaging.csv")

    print("Fetching DepMap model-level TPMs...")
    get_cache_item(["DepMap", "OmicsExpressionProteinCodingGenesTPMLogp1.csv"]).download(raw_dir / "OmicsExpressionProteinCodingGenesTPMLogp1.csv")

    print("Fetching DepMap gene-level copy number data...")
    get_cache_item(["DepMap", "OmicsCNGeneWGS.csv"]).download(raw_dir / "OmicsCNGeneWGS.csv")

    print("Fetching UniProt internal id mapping...")
    get_cache_item(["UniProt", "9606", "HUMAN_9606_idmapping.dat.gz"]).download(raw_dir / "HUMAN_9606_idmapping.dat.gz")
    uncompress(raw_dir / "HUMAN_9606_idmapping.dat.gz", raw_dir / "HUMAN_9606_idmapping.tsv")

    print("Fetching UniProt id external database mapping...")
    get_cache_item(["UniProt", "9606", "HUMAN_9606_idmapping_selected.tab.gz"]).download(raw_dir / "HUMAN_9606_idmapping_selected.tab.gz")
    uncompress(raw_dir / "HUMAN_9606_idmapping_selected.tab.gz", raw_dir / "HUMAN_9606_idmapping_selected.tsv")

    print("Fetching UniProt SwissProt genes...")
    get_cache_item(["UniProt", "9606", "SwissProt_9606.tsv"]).download(raw_dir / "SwissProt_9606.tsv")


if __name__ == "__main__":
    main()
