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
(dir_path / 'logs').mkdir(exist_ok=True)
logger.add(dir_path / 'logs' / "cache.log")

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
            with open(output, 'wb') as f:
                shutil.copyfileobj(response.raw, f)
            return response
    
    @staticmethod
    def coerce(obj: 'Service | str') -> 'Service':
        # TODO: This could also be replaced by coercing str to Service in CacheItem via pydantic.
        if isinstance(obj, str): return Service(url=obj)
        else: return obj

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

        with NamedTemporaryFile() as cached_file:
            logger.info(f"Downloading cache {self.cached}...")
            gdown.download(self.cached, cached_file)

            if self.pinned is not None:
                logger.info(f"Downloading pinned URL {self.pinned}...")
                Service.coerce(self.pinned).download(output)

                logger.info("Checking that the downloaded pinned artifact matches with cached artifact...")
                assert filecmp.cmp(output, cached_file.name)
            
            if self.unpinned is not None:
                logger.info(f"Downloading unpinned URL {self.unpinned}...")
                with NamedTemporaryFile() as unpinned_file:
                    Service.coerce(self.unpinned).download(unpinned_file.name)

                    logger.info("Checking that the downloaded unpinned artifact matches with cached artifact...")
                    if not filecmp.cmp(unpinned_file.name, cached_file.name):
                        # This gets saved to a file. Search for `logger.add` for more info.
                        logger.warning(f"Unpinned file {self.unpinned} for {self.name} does not match cache - this source should be updated!")


CacheDirectory = dict[str, Union[CacheItem, "CacheDirectory"]]

# An *unversioned* directory list.
directory: CacheDirectory = {
    "STRING": {
        "9606": {
            "9606.protein.links.txt.gz": CacheItem(
                name="STRING 9606 protein links",
                cached="https://drive.google.com/uc?id=13tE_-A6g7McZs_lZGz9As7iE-5cBFvqE",
                pinned="http://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz",
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
    "iRefIndex": {
        # This can also be obtained from the SPRAS repo
        # (https://github.com/Reed-CompBio/spras/blob/b5d7a2499afa8eab14c60ce0f99fa7e8a23a2c64/input/phosphosite-irefindex13.0-uniprot.txt).
        # iRefIndex has been down for quite some time, so this is only from the cache.
        "phosphosite-irefindex13.0-uniprot.txt": CacheItem(
            name="iRefIndex v13.0 UniProt interactome",
            cached="https://drive.google.com/uc?id=1fQ8Z3FjEwUseEtsExO723zj7mAAtdomo"
        )
    },
    "OsmoticStress": {
        "yeast_pcsf_network.sif": CacheItem(
            # In the paper https://doi.org/10.1016/j.celrep.2018.08.085
            name="Case Study Edge Results, from Supplementary Data 3",
            cached="https://drive.google.com/uc?id=1Agte0Aezext-8jLhGP4GmaF3tS7gHX-h"
        ),
        # The following files are from https://github.com/gitter-lab/osmotic-stress.
        # While the following files do point to the repository's main branch,
        # they aren't expected to actually change.
        "prizes.txt": CacheItem(
            name="Osmotic Stress Prizes",
            pinned="https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/Input%20Data/prizes.txt",
            cached="https://drive.google.com/uc?id=16WDQs0Vjv6rI12-hbifsbnpH31jMGhJg"
        ),
        "ChasmanNetwork-DirUndir.txt": CacheItem(
            name="Network Input",
            pinned="https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/Input%20Data/ChasmanNetwork-DirUndir.txt",
            cached="https://drive.google.com/uc?id=1qYXPaWcPU72YYME7NaBzD7thYCHRzrLH"
        ),
        "dummy.txt": CacheItem(
            name="Dummy Nodes File",
            pinned="https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/Input%20Data/dummy.txt",
            cached="https://drive.google.com/uc?id=1dsFIhBrIEahggg0JPxw64JwS51pKxoQU"
        ),
        "_edgeFreq.eda ": CacheItem(
            name="Case Study Omics Integrator Edge Frequencies",
            pinned="https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/Notebooks/Forest-TPS/_edgeFreq.eda",
            cached="https://drive.google.com/uc?id=1M_rxEzUCo_EVuFyM47OEH2J-4LB3eeCR"
        ),
        "goldStandardUnionDetailed.txt": CacheItem(
            name="Gold Standard Reference Pathways",
            pinned="https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/data/evaluation/goldStandardUnionDetailed.txt",
            cached="https://drive.google.com/uc?id=1-_zF9oKFCNmJbDCC2vq8OM17HJw80s2T"
        ),
    },
    "EGFR": {
        # The following files are from https://github.com/gitter-lab/tps.
        # While the following files do point to the repository's main branch,
        # they aren't expected to actually change.
        "eight-egfr-reference-all.txt": CacheItem(
            name="EGFR Gold Standard Reference",
            pinned="https://raw.githubusercontent.com/gitter-lab/tps/refs/heads/master/data/resources/eight-egfr-reference-all.txt",
            cached="https://drive.google.com/uc?id=15MqpIbH1GRA1tq0ZXH9oMnKytoFSzXyw"
        ),
        "egfr-prizes.txt": CacheItem(
            name="EGFR prizes",
            pinned="https://raw.githubusercontent.com/gitter-lab/tps/refs/heads/master/data/pcsf/egfr-prizes.txt",
            cached="https://drive.google.com/file/d/1nI5hw-rYRZPs15UJiqokHpHEAabRq6Xj/view?usp=sharing"
        )
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

    return current_item
