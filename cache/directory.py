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


@dataclass
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

        # TODO: yikes! same with self.unpinned


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
    "DISEASES": {
        # Instead of going through https://unmtid-shinyapps.net/shiny/tiga/, we use their
        # archived files directory instead.
        "tiga_gene-trait_stats.tsv": CacheItem(
            name="TIGA data",
            cached="https://drive.google.com/uc?id=114qyuNDy4qdmYDHHJAW-yBeTxcGTDUnK",
            pinned="https://unmtid-dbs.net/download/TIGA/20250916/tiga_gene-trait_stats.tsv",
        ),
        "HumanDO.tsv": CacheItem(
            name="Disease ontology data",
            cached="https://drive.google.com/uc?id=1lfB1DGJgrXTxP_50L6gGu_Nq6OyDjiIi",
            # DiseaseOntology is a decently updating repository!
            unpinned="https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/refs/heads/main/DOreports/HumanDO.tsv",
        ),
        "human_disease_textmining_filtered.tsv": CacheItem(
            name="DISEASES textmining channel",
            cached="https://drive.google.com/uc?id=1vD8KbT9sk04VEJx9r3_LglCTGYJdhN0D",
            unpinned="https://download.jensenlab.org/human_disease_textmining_filtered.tsv",
        ),
        "human_disease_knowledge_filtered.tsv": CacheItem(
            name="DISEASES knowledge channel",
            cached="https://drive.google.com/uc?id=1qGUnjVwF9-8p5xvp8_6CfVsbMSM_wkld",
            unpinned="https://download.jensenlab.org/human_disease_knowledge_filtered.tsv",
        ),
    },
    "BioMart": {
        "ensg-ensp.tsv": CacheItem(
            name="BioMart ENSG <-> ENSP mapping",
            cached="https://drive.google.com/uc?id=1-gPrDoluXIGydzWKjWEnW-nWhYu3YkHL",
            unpinned=fetch_biomart_service((dir_path / "biomart" / "ensg-ensp.xml").read_text()),
        )
    },
    "DepMap": {
        "OmicsProfiles.csv": CacheItem(
            name="DepMap omics metadata",
            cached="https://drive.google.com/uc?id=1i54aKfO0Ci2QKLTNJnuQ_jgGhH4c9rTL",
            pinned="https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2F2025-05-01-master-mapping-table-28c2.12%2Fpublic_release_date.2025-05-01.master_mapping_table.csv&dl_name=OmicsProfiles.csv&bucket=depmap-external-downloads",
        ),
        "CRISPRGeneDependency.csv": CacheItem(
            name="DepMap gene dependency probability estimates",
            cached="https://drive.google.com/uc?id=122rWNqT_u3M7B_11WYZMtOLiPbBykkaz",
            pinned="https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2F25q2-public-557c.3%2FCRISPRGeneDependency.csv&dl_name=CRISPRGeneDependency.csv&bucket=depmap-external-downloads",
        ),
        "OmicsSomaticMutationsMatrixDamaging.csv": CacheItem(
            name="DepMap genotyped matrix",
            cached="https://drive.google.com/uc?id=1W7N2H0Qi7NwmTmNChcwa2ZZ4WxAuz-Xh",
            pinned="https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2Fpublic-25q2-c5ef.87%2FOmicsSomaticMutationsMatrixDamaging.csv&dl_name=OmicsSomaticMutationsMatrixDamaging.csv&bucket=depmap-external-downloads",
        ),
        "OmicsExpressionProteinCodingGenesTPMLogp1.csv": CacheItem(
            name="DepMap model-level TPMs",
            cached="https://drive.google.com/uc?id=1P0m88eXJ8GPdru8h9oOcHPeXKU7ljIrP",
            pinned="https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2Fpublic-25q2-c5ef.73%2FOmicsExpressionProteinCodingGenesTPMLogp1.csv&dl_name=OmicsExpressionProteinCodingGenesTPMLogp1.csv&bucket=depmap-external-downloads",
        ),
        "OmicsCNGeneWGS.csv": CacheItem(
            name="DepMap gene-level copy number data",
            cached="https://drive.google.com/uc?id=1TPp3cfK7OZUrftucr3fLO-krXSQAA6Ub",
            pinned="https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2Fpublic-25q2-c5ef.104%2FOmicsCNGeneWGS.csv&dl_name=OmicsCNGeneWGS.csv&bucket=depmap-external-downloads",
        ),
    },
    "KEGG": {
        # For some reason, KEGG requires a Referer header: opening this URL otherwise fails.
        "ko03250.xml": CacheItem(
            name="KEGG 03250",
            cached="https://drive.google.com/uc?id=16dtWKHCQMp2qrLfFDE7nVhbwBCr2H5a9",
            unpinned=Service(
                "https://www.kegg.jp/kegg-bin/download?entry=ko03250&format=kgml", headers={"Referer": "https://www.kegg.jp/pathway/ko03250"}
            ),
        )
    },
    "HIV1": {
        # The following files are from https://github.com/gitter-lab/hiv1-aurkb.
        # While the following files do point to the repository's main branch,
        # they aren't expected to actually change.
        "prize_05.tsv": CacheItem(
            name="HIV_05 prizes",
            cached="https://drive.google.com/uc?id=1jVWNRPfYkbqimO44GdzXYB3-7NXhet1m",
            pinned="https://raw.githubusercontent.com/gitter-lab/hiv1-aurkb/refs/heads/main/Results/base_analysis/prize_05.csv",
        ),
        "prize_060.tsv": CacheItem(
            name="HIV_060 prizes",
            cached="https://drive.google.com/uc?id=1Aucgp7pcooGr9oT4m2bvYEuYW6186WxQ",
            pinned="https://raw.githubusercontent.com/gitter-lab/hiv1-aurkb/refs/heads/main/Results/base_analysis/prize_060.csv",
        ),
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
    "OsmoticStress": {
        "yeast_pcsf_network.sif": CacheItem(
            # In the paper https://doi.org/10.1016/j.celrep.2018.08.085
            name="Case Study Edge Results, from Supplementary Data 3",
            cached="https://drive.google.com/uc?id=1Agte0Aezext-8jLhGP4GmaF3tS7gHX-h",
        ),
        # The following files are from https://github.com/gitter-lab/osmotic-stress.
        # While the following files do point to the repository's main branch,
        # they aren't expected to actually change.
        "prizes.txt": CacheItem(
            name="Osmotic Stress Prizes",
            pinned="https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/Input%20Data/prizes.txt",
            cached="https://drive.google.com/uc?id=16WDQs0Vjv6rI12-hbifsbnpH31jMGhJg",
        ),
        "ChasmanNetwork-DirUndir.txt": CacheItem(
            name="Network Input",
            pinned="https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/Input%20Data/ChasmanNetwork-DirUndir.txt",
            cached="https://drive.google.com/uc?id=1qYXPaWcPU72YYME7NaBzD7thYCHRzrLH",
        ),
        "dummy.txt": CacheItem(
            name="Dummy Nodes File",
            pinned="https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/Input%20Data/dummy.txt",
            cached="https://drive.google.com/uc?id=1dsFIhBrIEahggg0JPxw64JwS51pKxoQU",
        ),
        "_edgeFreq.eda ": CacheItem(
            name="Case Study Omics Integrator Edge Frequencies",
            pinned="https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/Notebooks/Forest-TPS/_edgeFreq.eda",
            cached="https://drive.google.com/uc?id=1M_rxEzUCo_EVuFyM47OEH2J-4LB3eeCR",
        ),
        "goldStandardUnionDetailed.txt": CacheItem(
            name="Gold Standard Reference Pathways",
            pinned="https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/data/evaluation/goldStandardUnionDetailed.txt",
            cached="https://drive.google.com/uc?id=1-_zF9oKFCNmJbDCC2vq8OM17HJw80s2T",
        ),
    },
    "EGFR": {
        # The following files are from https://github.com/gitter-lab/tps.
        # While the following files do point to the repository's main branch,
        # they aren't expected to actually change.
        "eight-egfr-reference-all.txt": CacheItem(
            name="EGFR Gold Standard Reference",
            pinned="https://raw.githubusercontent.com/gitter-lab/tps/refs/heads/master/data/resources/eight-egfr-reference-all.txt",
            cached="https://drive.google.com/uc?id=15MqpIbH1GRA1tq0ZXH9oMnKytoFSzXyw",
        ),
        "egfr-prizes.txt": CacheItem(
            name="EGFR prizes",
            pinned="https://raw.githubusercontent.com/gitter-lab/tps/refs/heads/master/data/pcsf/egfr-prizes.txt",
            cached="https://drive.google.com/uc?id=1nI5hw-rYRZPs15UJiqokHpHEAabRq6Xj",
        ),
    },
    "Surfaceome": {
        "table_S3_surfaceome.xlsx": CacheItem(
            name="Human surfaceome",
            unpinned="http://wlab.ethz.ch/surfaceome/table_S3_surfaceome.xlsx",
            cached="https://docs.google.com/uc?id=1cBXYbDnAJVet0lv3BRrizV5FuqfMbBr0",
        )
    },
    "TranscriptionFactors": {
        "Homo_sapiens_TF.tsv": CacheItem(
            name="Human transcription factors",
            # This server has anti-bot protection, so to respect their wishes, we don't download from the server.
            # The original URL is https://guolab.wchscu.cn/AnimalTFDB4_static/download/TF_list_final/Homo_sapiens_TF,
            # which is accessible from https://guolab.wchscu.cn/AnimalTFDB4//#/Download -> Homo sapiens
            # (also under the Internet Archive as of Feb 2nd, 2026. If the original artifact disappears, the drive link below should suffice.)
            cached="https://drive.google.com/uc?id=1fVi18GpudUlquRPHgUJl3H1jy54gO-uz",
        )
    },
    "PathwayCommons": {
        "pc-biopax.owl.gz": CacheItem(
            name="PathwayCommons Universal BioPAX file",
            cached="https://drive.google.com/uc?id=1R7uE2ky7fGlZThIWCOblu7iqbpC-aRr0",
            pinned="https://download.baderlab.org/PathwayCommons/PC2/v14/pc-biopax.owl.gz",
        ),
        "pathways.txt.gz": CacheItem(
            name="PathwayCommons Pathway Identifiers",
            cached="https://drive.google.com/uc?id=1SMwuuohuZuNFnTev4zRNJrBnBsLlCHcK",
            pinned="https://download.baderlab.org/PathwayCommons/PC2/v14/pathways.txt.gz",
        ),
        "denylist.txt": CacheItem(
            name="PathwayCommons small molecule denylist",
            cached="https://drive.google.com/uc?id=1QmISJXPvVljA8oKuNYRUNbJJvZKPa_-u",
            pinned="https://download.baderlab.org/PathwayCommons/PC2/v14/blacklist.txt",
        ),
        "intermediate": {
            "pc-panther-biopax.owl": CacheItem(
                name="PathwayCommons PANTHER-only BioPAX file", cached="https://drive.google.com/uc?id=1MklrD8CJ1BIjh_wWr_g5rrIJ5XJB7FUI"
            )
        },
    },
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

    # Google Drive validation. TODO: remove if move to OSDF.
    if "uc?id=" not in current_item.cached or "/view?usp=sharing" in current_item.cached:
        raise RuntimeError(
            "Make sure your Google Drive URLs are in https://drive.google.com/uc?id=... format "
            + "with no /view?usp=sharing at the end. See CONTRIBUTING.md for more info."
        )

    return current_item
