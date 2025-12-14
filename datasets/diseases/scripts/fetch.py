"""
Fetches the latest DISEASES database channels, TIGA data, and human disease ontology data that we need.

Download pages:
- DISEASES: https://diseases.jensenlab.org/Downloads
- TIGA: https://unmtid-shinyapps.net/shiny/tiga/
- Disease Ontology: https://disease-ontology.org/downloads/
"""

from pathlib import Path
import os
import urllib.request

# https://stackoverflow.com/a/5137509/7589775
dir_path = os.path.dirname(os.path.realpath(__file__))

raw_dir = Path(dir_path, "..", "raw")

def main():
    # We only need the text mining and knowledge channels
    # and avoid the integrated channel as it is the multiplied probabilities of all
    # three channels (as mentioned from Damian Szklarczyk)

    raw_dir.mkdir(exist_ok=True)

    print("Fetching DISEASES text channel...")
    urllib.request.urlretrieve(
        "https://download.jensenlab.org/human_disease_textmining_filtered.tsv",
        raw_dir / "human_disease_textmining_filtered.tsv"
    )

    print("Fetching DISEASES knowledge channel...")
    urllib.request.urlretrieve(
        "https://download.jensenlab.org/human_disease_knowledge_filtered.tsv",
        raw_dir / "human_disease_knowledge_filtered.tsv"
    )

    print("Fetching TIGA data...")
    urllib.request.urlretrieve(
        "https://unmtid-shinyapps.net/shiny/tiga/session/4fc1d570b85bb6fcc4e9660a7944a6e3/download/gt_file?w=",
        raw_dir / "tiga_gene-trait_stats.tsv"
    )

    print("Fetching human disease ontology data...")
    urllib.request.urlretrieve(
        "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/016a4ec33d1a1508d669650086cd92ccebe138e6/DOreports/HumanDO.tsv",
        raw_dir / "HumanDO.tsv"
    )

    # TODO: we should make some assertions about these file structures to check for correctness early

if __name__ == '__main__':
    main()
