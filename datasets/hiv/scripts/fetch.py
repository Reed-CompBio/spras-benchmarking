"""
Fetches the `prizes_05.tsv` and `prizes_060.tsv` files from
https://github.com/gitter-lab/hiv1-aurkb, as well as the KGML file
for `kegg_orthology.py`

Associated paper: https://doi.org/10.1371/journal.ppat.1011492
"""

import urllib.request
import os
from pathlib import Path

# https://stackoverflow.com/a/5137509/7589775
hiv_path = Path(os.path.dirname(os.path.realpath(__file__))).parent

base_url = "https://raw.githubusercontent.com/gitter-lab/hiv1-aurkb/ac9278d447e4188eea3bf4b24c4c4e0c19b0c6d9/Results/base_analysis/"
prizes_05_url = base_url + "prize_05.csv"
prizes_060_url = base_url + "prize_060.csv"

def main():
    # Note: These files are .tsv, but have the wrong file extension .csv in the original data source.
    urllib.request.urlretrieve(prizes_05_url, hiv_path / "raw" / "prizes_05.tsv")
    urllib.request.urlretrieve(prizes_060_url, hiv_path / "raw" / "prizes_060.tsv")

    # and our final KGML file for the HIV pathway.
    # KEGG requires a Referrer (server-CORS enforcement?)
    # https://stackoverflow.com/a/46511429/7589775
    opener = urllib.request.build_opener()
    opener.addheaders = [('Referer', 'https://www.kegg.jp/pathway/ko03250')]
    urllib.request.install_opener(opener)
    urllib.request.urlretrieve(
        "https://www.kegg.jp/kegg-bin/download?entry=ko03250&format=kgml",
        hiv_path / "raw" / "ko03250.xml")

if __name__ == '__main__':
    main()
