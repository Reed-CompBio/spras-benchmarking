from dataclasses import dataclass
from typing import Union
from os import PathLike
from tempfile import NamedTemporaryFile
import urllib.request
import filecmp

import gdown

@dataclass
class CacheItem:
    """Class for differentriating between offline and online items in a cache."""
    cached: str
    online: str

    def download(self, output: str | PathLike):
        print(f"Downloading {self.online}...")

        urllib.request.urlretrieve(self.online, output)

        with NamedTemporaryFile() as cached_file:
            print(f"Downloading cache {self.cached}...")
            gdown.download(self.cached, cached_file)
            print("Checking that downloaded artifact matches with cached artifact...")
            filecmp.cmp(output, cached_file.name)

CacheDirectory = dict[str, Union[CacheItem, 'CacheDirectory']]

# An *unversioned* directory list.
directory: CacheDirectory = {
    'STRING': {
        '9606': CacheItem(
            cached='https://drive.google.com/uc?id=1fvjdIbgzbgJrdJxWRRRwwS1zuegf6DOj',
            online='http://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz'
        )
    },
    'DISEASES': {
        'tiga_gene-trait_stats.tsv': CacheItem(
            cached='https://drive.google.com/uc?id=114qyuNDy4qdmYDHHJAW-yBeTxcGTDUnK',
            online='https://unmtid-shinyapps.net/shiny/tiga/session/3540ae3a93d3494ca59ebe1ace3fc6be/download/gt_file?w='
        ),
        'HumanDO.tsv': CacheItem(
            cached='https://drive.google.com/uc?id=1lfB1DGJgrXTxP_50L6gGu_Nq6OyDjiIi',
            online='https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/016a4ec33d1a1508d669650086cd92ccebe138e6/DOreports/HumanDO.tsv'
        ),
        'human_disease_textmining_filtered.tsv': CacheItem(
            cached='https://drive.google.com/uc?id=1vD8KbT9sk04VEJx9r3_LglCTGYJdhN0D',
            online='https://download.jensenlab.org/human_disease_textmining_filtered.tsv'
        ),
        'human_disease_knowledge_filtered.tsv': CacheItem(
            cached='https://drive.google.com/uc?id=1qGUnjVwF9-8p5xvp8_6CfVsbMSM_wkld',
            online='https://download.jensenlab.org/human_disease_knowledge_filtered.tsv'
        ),
    }
}

def get_cache_item(path: list[str]) -> CacheItem:
    """Takes a path and gets the underlying cache item."""
    assert len(path) != 0

    current_item = directory
    for entry in path:
        if isinstance(current_item, CacheItem):
            raise ValueError(f"Path {path} leads to a cache item too early!")
        current_item = current_item[entry]
    
    if not isinstance(current_item, CacheItem):
        raise ValueError(f"Path {path} doesn't lead to a cache item")

    return current_item
