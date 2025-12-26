"""
Fetches the latest DISEASES database channels, TIGA data, and human disease ontology data that we need.

Download pages:
- DISEASES: https://diseases.jensenlab.org/Downloads
- TIGA: https://unmtid-shinyapps.net/shiny/tiga/
- Disease Ontology: https://disease-ontology.org/downloads/
"""

from pathlib import Path
import os
from cache.directory import get_cache_item

# https://stackoverflow.com/a/5137509/7589775
dir_path = os.path.dirname(os.path.realpath(__file__))

raw_dir = Path(dir_path, "..", "raw")


def main():
    # We only need the text mining and knowledge channels
    # and avoid the integrated channel as it is the multiplied probabilities of all
    # three channels (personal correspondence with Damian Szklarczyk)

    raw_dir.mkdir(exist_ok=True)

    print("Fetching DISEASES text channel...")
    get_cache_item(["DISEASES", "human_disease_textmining_filtered.tsv"]).download(raw_dir / "human_disease_textmining_filtered.tsv")

    print("Fetching DISEASES knowledge channel...")
    get_cache_item(["DISEASES", "human_disease_knowledge_filtered.tsv"]).download(raw_dir / "human_disease_knowledge_filtered.tsv")

    print("Fetching TIGA data...")
    get_cache_item(["DISEASES", "tiga_gene-trait_stats.tsv"]).download(raw_dir / "tiga_gene-trait_stats.tsv")

    print("Fetching human disease ontology data...")
    get_cache_item(["DISEASES", "HumanDO.tsv"]).download(raw_dir / "HumanDO.tsv")

    print("Fetching BioMart ENSG - ENSP mapping...")
    get_cache_item(["BioMart", "ensg-ensp.tsv"]).download(raw_dir / "ensg-ensp.tsv")


if __name__ == "__main__":
    main()
