"""Fetches the generic files needed for the WikiPathways dataset."""

from pathlib import Path
import os
from cache.directory import get_cache_item

dir_path = os.path.dirname(os.path.realpath(__file__))

raw_dir = Path(dir_path, "..", "raw")


def main():
    raw_dir.mkdir(exist_ok=True)

    print("Fetching BioMart Entrez -> ENSG mapping...")
    get_cache_item(["BioMart", "entrez-ensg.tsv"]).download(raw_dir / "entrez-ensg.tsv")

    # See special_nodes.py on how these files are taken from the respective supplementary info sections.
    print("Fetching supplementary info from https://doi.org/10.1186/1741-7007-7-50:")
    get_cache_item(["DOI", "10.1186/1741-7007-7-50", "12915_2009_258_MOESM1_ESM - Data.tsv"]).download(raw_dir / "10_1186-data.tsv")

    print("Fetching supplementary info from https://doi.org/10.1016/j.cell.2010.01.044:")
    get_cache_item(["DOI", "10.1016/j.cell.2010.01.044", "mmc1.tsv"]).download(raw_dir / "10_1016-mmc1.tsv")

    print("Fetching supplementary info from https://doi.org/10.1038/nrg2538:")
    get_cache_item(["DOI", "10.1038/nrg2538", "st2.tsv"]).download(raw_dir / "st2.tsv")


if __name__ == "__main__":
    main()
