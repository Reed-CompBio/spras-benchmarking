"""
Fetches the second supplementary data from the TPS paper:
https://doi.org/10.1016/j.celrep.2018.08.085. We trust the processed data from the paper:
see the README for motivation.

The ZIP contains `p-values-first.tsv` and `p-values-prev.tsv`,
which are fed in to generate the prizes file.
"""

import requests
import io
from pathlib import Path
import os
import zipfile

current_directory = Path(os.path.dirname(os.path.realpath(__file__)))

DOWNLOAD_URL = "https://ars.els-cdn.com/content/image/1-s2.0-S2211124718313895-mmc3.zip"

def main():
    # https://stackoverflow.com/a/14260592/7589775
    request = requests.get(DOWNLOAD_URL)
    zipf = zipfile.ZipFile(io.BytesIO(request.content))

    download_folder = current_directory / '..' / 'download'
    download_folder.mkdir(exist_ok=True)
    zipf.extractall(current_directory / '..' / 'download')

if __name__ == '__main__':
    main()
