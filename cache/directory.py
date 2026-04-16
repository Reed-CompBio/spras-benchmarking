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

# Our cache emits warnings for files with unpinned versions that don't match the cache.
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
        with requests.get(self.url, stream=True, headers=self.headers) as response:
            response.raw.decode_content = True
            with open(output, "wb") as f:
                shutil.copyfileobj(response.raw, f)
            return response

    # NOTE: this is slightly yucky code deduplication. The only intended values of `downloaded_file_type` are `pinned` and `unpinned`.
    def download_against_cache(self, cache: Path, downloaded_file_type: str, move_output: bool):
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
            if move_output:
                shutil.move(cache, debug_file_path)
            else:
                shutil.copy(cache, debug_file_path)
            # We use a custom error type to prevent any overlap with RuntimeError. I am not sure if there is any.
            raise DownloadFileCheckException(
                f"The {downloaded_file_type} file {downloaded_file_path} and "
                + f"cached file originally at {cache} do not match! "
                + f"Compare the pinned {downloaded_file_path} and the cached {debug_file_path}."
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


def fetch_biomart_service(xml: str) -> Service:
    """
    Access BioMart data through the BioMart REST API:
    https://useast.ensembl.org/info/data/biomart/biomart_restful.html#biomartxml
    """
    ROOT = "http://www.ensembl.org/biomart/martservice?query="
    return Service(ROOT + urllib.parse.quote_plus(xml))


@dataclass(frozen=True)
class CacheItem:
    """
    Class for differentriating between different ways of fetching data.
    As mentioned in the ./README.md, `cached` is always needed, and we differentriate between service outage (`pinned`)
    and data needing updates (`unpinned`). There is no need to specify both keys at once, but the choice does matter
    for how errors are presented during benchmarking runs.
    """

    name: str
    """The display name of the artifact, used for human-printing."""

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
    we say that the file has a new version.

    If `pinned` is None and `unpinned` doesn't match `cached`, we warn instead of erroring.

    We will still error if the status code is not 2XX (a successful request).
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
        logger.info(f"Fetching {self.name}...")

        logger.info(f"Downloading cache {self.cached} to {output}...")
        gdown.download(self.cached, str(output))  # gdown doesn't have a type signature, but it expects a string :/

        if self.pinned is not None:
            Service.coerce(self.pinned).download_against_cache(cache=Path(output), downloaded_file_type="pinned", move_output=True)
        if self.unpinned is not None:
            # Normally, download_against_cache raises a DownloadFileCheckException: we catch it and warn instead if that happens.
            try:
                Service.coerce(self.unpinned).download_against_cache(cache=Path(output), downloaded_file_type="unpinned", move_output=False)
            except DownloadFileCheckException as err:
                logger.warning(err)


CacheDirectory = dict[str, Union[CacheItem, "CacheDirectory"]]

# An *unversioned* directory list.
directory: CacheDirectory = {
    "STRING": {
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
    "UniProt": {
        # We use FTP when possible, but we delegate to the UniProt REST API in cases that would save significant bandwidth.
        # See https://ftp.uniprot.org/pub/databases/uniprot/current_release/README for the FTP README.
        "9606": {
            # We prefer manually curated, or SwissProt, genes. This URL selects these genes using the REST API.
            "SwissProt_9606.tsv": CacheItem(
                name="UniProt 9606 SwissProt genes",
                cached="https://drive.google.com/uc?id=1h2Cl-60qcKse-djcsqlRXm_n60mVY7lk",
                unpinned="https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names&format=tsv&query=%28*%29+AND+%28reviewed%3Atrue%29+AND+%28model_organism%3A9606%29",
            ),
            # Sources
            "sources.tsv": CacheItem(
                # Where KW-0675 is the UniProt keyword for receptors
                name="UniProt-tagged sources (receptors)",
                unpinned="https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid&format=tsv&query=%28%28keyword%3A%22KW-0675%22%29%29+AND+%28reviewed%3Atrue%29+AND+%28model_organism%3A9606%29",
                cached="https://drive.google.com/uc?id=1VbCLH9yoJ41QhzhsSy9ICAU2MLAAxfJe"
            ),
            # Targets
            "targets.tsv": CacheItem(
                name="UniProt-tagged targets (transcription factors)",
                # Where KW-0539 and KW-0805 are the UniProt keywords for the nucleus and transcription regulators, respectively.
                unpinned="https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid&format=tsv&query=%28%28keyword%3AKW-0539%29+OR+%28keyword%3AKW-0805%29%29+AND+%28reviewed%3Atrue%29+AND+%28model_organism%3A9606%29",
                cached="https://drive.google.com/uc?id=1gg_2IO1xHeho8KkcYVIfqHNWSRZx6gd1"
            ),
            # idmapping FTP files. See the associated README:
            # https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
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
    "BioMart": {
        "ensg-ensp.tsv": CacheItem(
            name="BioMart ENSG <-> ENSP mapping",
            cached="https://drive.google.com/uc?id=1-gPrDoluXIGydzWKjWEnW-nWhYu3YkHL",
            unpinned=fetch_biomart_service((dir_path / "biomart" / "ensg-ensp.xml").read_text()),
        )
    },
    "iRefIndex": {
        # This can also be obtained from the SPRAS repo, though the SPRAS repo removes self loops. We don't.
        # (https://github.com/Reed-CompBio/spras/blob/b5d7a2499afa8eab14c60ce0f99fa7e8a23a2c64/input/phosphosite-irefindex13.0-uniprot.txt).
        # iRefIndex has been down for quite some time, so we grab this from a repository instead.
        # While the following files do point to the repository's main branch,
        # they aren't expected to actually change, so we make them `pinned`.
        "phosphosite-irefindex13.0-uniprot.txt": CacheItem(
            name="iRefIndex v13.0 UniProt interactome",
            cached="https://drive.google.com/uc?id=1fQ8Z3FjEwUseEtsExO723zj7mAAtdomo",
            pinned="https://raw.githubusercontent.com/gitter-lab/tps/refs/heads/master/data/networks/phosphosite-irefindex13.0-uniprot.txt",
        )
    },
    "PathwayCommons": {
        "pathways.txt.gz": CacheItem(
            name="PathwayCommons Pathway Identifiers",
            cached="https://drive.google.com/uc?id=1SMwuuohuZuNFnTev4zRNJrBnBsLlCHcK",
            pinned="https://download.baderlab.org/PathwayCommons/PC2/v14/pathways.txt.gz",
        ),
    },

    "DepMap_25Q3": {
        "OmicsProfiles_25Q3.csv": CacheItem(
            name="DepMap omics metadata",
            cached="https://drive.google.com/uc?id=10eJ6ZcDlxEBsfHsvkmrgtTEd__71fC_7",
            # pinned="https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2F2025-05-01-master-mapping-table-28c2.12%2Fpublic_release_date.2025-05-01.master_mapping_table.csv&dl_name=OmicsProfiles.csv&bucket=depmap-external-downloads",
        ),

        "CRISPRGeneDependency_25Q3.csv": CacheItem(
            name="DepMap crispr gene dependency probability estimates",
            cached="https://drive.google.com/uc?id=18hGSWlWEo_R6FmYziZKEISYCbDqrNty_",
            # pinned="https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2F25q2-public-557c.3%2FCRISPRGeneDependency.csv&dl_name=CRISPRGeneDependency.csv&bucket=depmap-external-downloads",
        ),

        "OmicsSomaticMutationsMatrixDamaging_25Q3.csv": CacheItem(
            name="DepMap genotyped matrix",
            cached="https://drive.google.com/uc?id=1M-ybvvKvGbdRhPLGjtBzKJviVecVXbDi",
            # pinned="https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2Fpublic-25q2-c5ef.87%2FOmicsSomaticMutationsMatrixDamaging_25Q3.csv&dl_name=OmicsSomaticMutationsMatrixDamaging_25Q3.csv&bucket=depmap-external-downloads",
        ),

        "Model_25Q3.csv": CacheItem(
            name="DepMap Metadata describing all cancer models/cell lines which are referenced by a dataset contained within the DepMap portal.",
            cached="https://drive.google.com/uc?id=1-26H3bvAAEt7UjlxpiK1cxORiAGJnxU9"
        ),
    },

    # TODO: update the ccle data once the new pr gets merged: https://github.com/cBioPortal/datahub/pull/2283
    "cbioportal": {
        "data_cna_cbioportal_ccle2019.txt": CacheItem(
            name="Copy number alterations from CCLE",
            cached="https://drive.google.com/uc?id=1C-OQu80Ptfy0-aBWg6nlULcQQvfzBBPp"
        ),
    },

    "CCLE_2019": {
        "Cell_lines_annotations_20181226_ccle2019.txt": CacheItem(
            name = "Cell line annotations",
            cached="https://drive.google.com/uc?id=1MDe_MizqSfaH58UMrwzhorXStMNCNNBy"
        )
    },

    "estimateNCA": {
        "consensus_tfa_march_6.tsv": CacheItem(
            name = "",
            cached="https://drive.google.com/uc?id=1haTqjyqkWoYTLw_6tubc0rJUiJIszC6v"
        ),

        "tfs_beyond_2sd_per_cell_line.csv": CacheItem(
            name="",
            cached="https://drive.google.com/uc?id=1vl6t_8bXfbyCYzaw9u2HpU430O_IMHWm"
        )
    }

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
