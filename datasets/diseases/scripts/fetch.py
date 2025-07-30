"""
Fetches the latest [desired] DISEASES database channels.

See the diseases download page: https://diseases.jensenlab.org/Downloads
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

    # TODO: we should make some assertions about its structure to check for correctness early

if __name__ == '__main__':
    main()
