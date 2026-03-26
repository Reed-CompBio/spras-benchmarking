from dataclasses import dataclass
from typing import Union
from os import PathLike
from tempfile import NamedTemporaryFile
from typing import Optional, Mapping
import filecmp
from pathlib import Path
import warnings
import requests
import shutil
import urllib.parse

import gdown
from loguru import logger

dir_path = Path(__file__).parent.resolve()

# Our cache emits warnings for files with unpinned versions that don't match the cache
# using loguru, and warnings are added to a local `logs` folder.
(dir_path / "logs").mkdir(exist_ok=True)
logger.add(dir_path / "logs" / "cache.log", level="WARNING")


class DownloadFileCheckException(RuntimeError):
    """See Service#download_against_cache for some motivation for this custom error"""


@dataclass
class Service:
    url: str
    headers: Optional[Mapping[str, str]] = None

    def download(self, output: str | PathLike) -> requests.Response:
        """
        Downloads a URL, returning the response (to be used with `with`) and modifying the output path.
        """
        # As per https://stackoverflow.com/a/39217788/7589775 to enable download streaming.
        with requests.get(self.url, stream=True, headers=self.headers, allow_redirects=True) as response:
            response.raw.decode_content = True
            with open(output, "wb") as f:
                shutil.copyfileobj(response.raw, f)
            return response

    # NOTE: this is slightly yucky code deduplication. The only intended values of `downloaded_file_type` are `pinned` and `unpinned`.
    def download_against_cache(self, cache: Path, downloaded_file_type: str, move_output_on_error: bool):
        """
        Downloads `this` Service and checks it against the provided `cache` at path. In logs,
        the file will be referred to as `downloaded_file_type`.

        @param move_output: Whether or not output should be irrecoverably moved instead of just copied.
        """
        logger.info(f"Downloading {downloaded_file_type} file {self.url} to check against with artifact at {cache}...")
        downloaded_file_path = Path(NamedTemporaryFile(delete=False).name)

        self.download(downloaded_file_path)
        logger.info(f"Checking that the {downloaded_file_type} artifact {downloaded_file_path} matches with cached artifact at {cache}...")

        if not filecmp.cmp(cache, downloaded_file_path):
            # This entire if-branch is debug shenanigans: we want to be able to easily compare our current cached file to the online file,
            # especially since some `Service`s have special errors that can make the request hard to compare in the browser.

            debug_file_path = Path(NamedTemporaryFile(prefix="spras-benchmarking-debug-artifact", delete=False).name)
            # We use shutil over Path#rename since temporary directories can be mounted to a different file system.
            if move_output_on_error:
                shutil.move(cache, debug_file_path)
            else:
                shutil.copy(cache, debug_file_path)
            # We use a custom error type to prevent any overlap with RuntimeError. I am not sure if there is any.
            raise DownloadFileCheckException(
                f"The {downloaded_file_type} file {downloaded_file_path} and "
                + f"cached file originally at {cache} do not match! "
                + f"Compare the pinned {downloaded_file_path} and the cached {debug_file_path}. "
                + "If this file updated, please update the underlying `cache` file to match."
            )
        else:
            # Since we don't clean up pinned_file_path for the above branch's debugging,
            # we need to clean it up here.
            downloaded_file_path.unlink()

    @staticmethod
    def coerce(obj: "Service | str") -> "Service":
        # TODO: This could also be replaced by coercing str to Service in CacheItem via pydantic.
        if isinstance(obj, str):
            return Service(url=obj)
        return obj


def fetch_biomart_service(xml: str, archived: bool=False) -> Service:
    """
    Access BioMart data through the BioMart REST API:
    https://useast.ensembl.org/info/data/biomart/biomart_restful.html#biomartxml

    We also provide links to the Ensembl archives for pinned files.
    """
    ROOT = "http://www.ensembl.org/biomart/martservice?query="
    ROOT_ARCHIVED = "http://sep2025.archive.ensembl.org/biomart/martservice?query="
    return Service((ROOT_ARCHIVED if archived else ROOT) + urllib.parse.quote_plus(xml))


@dataclass(frozen=True)
class CacheItem:
    """
    Class for differentriating between different ways of fetching data.
    As mentioned in the ./README.md, `cached` is always needed, and we differentiate between service outage (`pinned`)
    and data needing updates (`unpinned`). There is no need to specify both keys at once, but the choice does matter
    for how errors are presented during benchmarking runs.
    """

    name: str
    """The display name of the artifact, used for human-readable logs."""

    cached: str
    """
    The URL of the cached file, which is currently a Google Drive URL.
    """

    pinned: Optional[Service | str] = None
    """
    The Service (URL + headers) of the file, which is the 'pinned' file.
    By a pinned file, we say that the file has a dedicated version, and should not change.
    If this is None, we go for the `unpinned` file or `cached` if `unpinned` is None.
    """

    unpinned: Optional[Service | str] = None
    """
    Analogously to `pinned`, this is a Service (URL + headers) which is 'unpinned,'
    or lacks a dedicated version. When `pinned` matches `cached` but `unpinned` doesn't match `pinned`,
    we say that the file has a new version but won't be automatically updated to the new version.

    If unpinned` doesn't match `cached`, we emit a warning.
    """

    def __post_init__(self):
        # Google Drive validation. TODO: remove if move to OSDF.
        if "uc?id=" not in self.cached or "/view?usp=sharing" in self.cached:
            raise RuntimeError(
                "Make sure your Google Drive URLs are in https://drive.google.com/uc?id=... format "
                + "with no /view?usp=sharing at the end. See CONTRIBUTING.md for more info."
            )

    @classmethod
    @warnings.deprecated("Pending for removal after the CONTRIBUTING guide is updated.")
    def cache_only(cls, name: str, cached: str) -> "CacheItem":
        """Wrapper method to explicitly declare a CacheItem as cached only."""
        return cls(name=name, cached=cached)

    def download(self, output: str | PathLike):
        """
        Downloads this `CacheItem` to the desired `output`,
        comparing the `cached` file to the `pinned` and `unpinned` files,
        warning when `cached` doesn't match `unpinned`, and erroring when
        `cached` doesn't match `pinned`.

        The file from `cache` is the file that gets downloaded to `output`.
        """
        logger.info(f"Fetching {self.name}...")

        logger.info(f"Downloading cache {self.cached} to {output}...")
        gdown.download(self.cached, str(output))  # gdown doesn't have a type signature, but it expects a string

        # If the file is pinned, we move the file to make sure it doesn't accidentally get used in workflows,
        # and stop the entire workflow if something bad happens. The converse for unpinned is in the other branch.
        if self.pinned is not None:
            Service.coerce(self.pinned).download_against_cache(cache=Path(output), downloaded_file_type="pinned", move_output_on_error=True)
        if self.unpinned is not None:
            # Normally, download_against_cache raises a DownloadFileCheckException: we catch it and warn instead if that happens.
            try:
                Service.coerce(self.unpinned).download_against_cache(cache=Path(output), downloaded_file_type="unpinned", move_output_on_error=False)
            except DownloadFileCheckException as err:
                logger.warning(err)


CacheDirectory = dict[str, Union[CacheItem, "CacheDirectory"]]

# An *unversioned* directory list.
directory: CacheDirectory = {
    # STRINGDB: https://string-db.org/
    # You can see more information about these files at https://string-db.org/cgi/download.
    "STRING": {
        # 9606 is human.
        "9606": {
            "9606.protein.links.full.txt.gz": CacheItem(
                name="STRING 9606 full protein links",
                cached="https://drive.google.com/uc?id=13tE_-A6g7McZs_lZGz9As7iE-5cBFvqE",
                pinned="http://stringdb-downloads.org/download/protein.links.full.v12.0/9606.protein.links.full.v12.0.txt.gz",
            ),
            "9606.protein.aliases.txt.gz": CacheItem(
                name="STRING 9606 protein aliases",
                cached="https://drive.google.com/uc?id=1IWrQeTVCcw1A-jDk-4YiReWLnwP0S9bY",
                pinned="https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz",
            ),
        }
    },
    # https://www.uniprot.org/
    "UniProt": {
        # We use FTP when possible, but we delegate to the UniProt REST API in cases that would save significant bandwidth.
        # See https://ftp.uniprot.org/pub/databases/uniprot/current_release/README for the FTP README.
        # 9606 is human.
        "9606": {
            # We prefer manually curated, or SwissProt, genes.
            # This URL selects these genes using the REST API.
            # UniProt REST doesn't seem to have any way to version it, so we only provide the `unpinned` URL.
            "SwissProt_9606.tsv": CacheItem(
                name="UniProt 9606 SwissProt genes",
                cached="https://drive.google.com/uc?id=1h2Cl-60qcKse-djcsqlRXm_n60mVY7lk",
                unpinned="https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names&format=tsv&query=%28*%29+AND+%28reviewed%3Atrue%29+AND+%28model_organism%3A9606%29",
            ),
            # idmapping FTP files. See the associated README:
            # https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
            # We use these files as our primary source of identifier mapping.
            # Unfortunately, there are no accompanying `pinned` URLs for these, as the previous releases
            # contain files with data magnitudes higher than anything we process.
            "HUMAN_9606_idmapping_selected.tab.gz": CacheItem(
                name="UniProt 9606 ID external database mapping",
                cached="https://drive.google.com/uc?id=1Oysa5COq31H771rVeyrs-6KFhE3VJqoX",
                unpinned="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz",
            ),
            "HUMAN_9606_idmapping.dat.gz": CacheItem(
                name="UniProt 9606 internal id mapping",
                cached="https://drive.google.com/uc?id=1lGxrx_kGyNdupwIOUXzfIZScc7rQKP-O",
                unpinned="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
            ),
        }
    },
    # https://www.ensembl.org/info/data/biomart/index.html
    "BioMart": {
        # ENSG to ENSP mapping
        # we generally prefer UniProt mappings over this, since we usually work over UniProt-compatible
        # ENSG/ENSP identifiers, but sometimes identifiers live outside of this space.
        "ensg-ensp.tsv": CacheItem(
            name="BioMart ENSG <-> ENSP mapping",
            cached="https://drive.google.com/uc?id=1-gPrDoluXIGydzWKjWEnW-nWhYu3YkHL",
            unpinned=fetch_biomart_service((dir_path / "biomart" / "ensg-ensp.xml").read_text()),
            pinned=fetch_biomart_service((dir_path / "biomart" / "ensg-ensp.xml").read_text(), archived=True),
        )
    },
    # https://www.pathwaycommons.org/
    "PathwayCommons": {
        "pathways.txt.gz": CacheItem(
            name="PathwayCommons Pathway Identifiers",
            cached="https://drive.google.com/uc?id=1SMwuuohuZuNFnTev4zRNJrBnBsLlCHcK",
            pinned="https://download.baderlab.org/PathwayCommons/PC2/v14/pathways.txt.gz",
        ),
    },
}


def get_cache_item(path: tuple[str, ...]) -> CacheItem:
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
