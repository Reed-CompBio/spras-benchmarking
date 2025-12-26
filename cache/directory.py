from dataclasses import dataclass
from typing import Union
from os import PathLike
from tempfile import NamedTemporaryFile
import urllib.request
import filecmp
import urllib.parse
import os
from pathlib import Path

import gdown

dir_path = Path(os.path.dirname(os.path.realpath(__file__)))

def fetch_biomart_url(xml: str) -> str:
    """
    Access BioMart data through the BioMart REST API:
    https://useast.ensembl.org/info/data/biomart/biomart_restful.html#biomartxml
    """
    ROOT = "http://www.ensembl.org/biomart/martservice?query="
    return ROOT + urllib.parse.quote_plus(xml)

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


CacheDirectory = dict[str, Union[CacheItem, "CacheDirectory"]]

# An *unversioned* directory list.
directory: CacheDirectory = {
    "STRING": {
        "9606": {
            "links": CacheItem(
                cached="https://drive.google.com/uc?id=1fvjdIbgzbgJrdJxWRRRwwS1zuegf6DOj",
                online="http://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz",
            ),
            "aliases": CacheItem(
                cached="https://drive.google.com/uc?id=1IWrQeTVCcw1A-jDk-4YiReWLnwP0S9bY",
                online="https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz",
            )
        }
    },
    "DISEASES": {
        # Instead of going through https://unmtid-shinyapps.net/shiny/tiga/, we use their
        # archived files directory instead.
        "tiga_gene-trait_stats.tsv": CacheItem(
            cached="https://drive.google.com/uc?id=114qyuNDy4qdmYDHHJAW-yBeTxcGTDUnK",
            online="https://unmtid-dbs.net/download/TIGA/20250916/tiga_gene-trait_stats.tsv",
        ),
        "HumanDO.tsv": CacheItem(
            cached="https://drive.google.com/uc?id=1lfB1DGJgrXTxP_50L6gGu_Nq6OyDjiIi",
            online="https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/016a4ec33d1a1508d669650086cd92ccebe138e6/DOreports/HumanDO.tsv",
        ),
        "human_disease_textmining_filtered.tsv": CacheItem(
            cached="https://drive.google.com/uc?id=1vD8KbT9sk04VEJx9r3_LglCTGYJdhN0D",
            online="https://download.jensenlab.org/human_disease_textmining_filtered.tsv",
        ),
        "human_disease_knowledge_filtered.tsv": CacheItem(
            cached="https://drive.google.com/uc?id=1qGUnjVwF9-8p5xvp8_6CfVsbMSM_wkld",
            online="https://download.jensenlab.org/human_disease_knowledge_filtered.tsv",
        ),
    },
    "BioMart": {
        "ensg-ensp.tsv": CacheItem(
            cached="https://drive.google.com/uc?id=1-gPrDoluXIGydzWKjWEnW-nWhYu3YkHL",
            online=fetch_biomart_url((dir_path / "biomart" / "ensg-ensp.xml").read_text())
        )
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
