"""Fetches the generic files needed for the WikiPathways dataset."""

from pathlib import Path
import os
from cache.directory import get_cache_item

dir_path = os.path.dirname(os.path.realpath(__file__))

raw_dir = Path(dir_path, "..", "raw")


def main():
    raw_dir.mkdir(exist_ok=True)

    print(f"Fetching BioMart Entrez -> ENSG mapping...")
    get_cache_item(["BioMart", "entrez-ensg.tsv"]).download(raw_dir / "entrez-ensg.tsv")

if __name__ == '__main__':
    main()
