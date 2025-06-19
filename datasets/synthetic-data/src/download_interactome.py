from pathlib import Path
import os
import urllib.request 

current_directory = Path(os.path.dirname(os.path.realpath(__file__)))
DOWNLOAD_LINK = "https://stringdb-downloads.org/download/protein.links.full.v12.0/9606.protein.links.full.v12.0.txt.gz"

if __name__ == '__main__':
    urllib.request.urlretrieve(DOWNLOAD_LINK, current_directory / '..' / 'intermediate' / '9606.protein.links.full.v12.0.txt.gz')
