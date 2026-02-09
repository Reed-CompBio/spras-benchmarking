from dataclasses import dataclass
from typing import Union
from os import PathLike
from tempfile import NamedTemporaryFile
from typing import Optional, Mapping
import filecmp
from pathlib import Path
from enum import Enum
import warnings
import requests
import shutil
import urllib.parse

import gdown

dir_path = Path(__file__).parent.resolve()

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

def fetch_biomart_service(xml: str) -> Service:
    """
    Access BioMart data through the BioMart REST API:
    https://useast.ensembl.org/info/data/biomart/biomart_restful.html#biomartxml
    """
    ROOT = "http://www.ensembl.org/biomart/martservice?query="
    return Service(ROOT + urllib.parse.quote_plus(xml))

class OnlineStatus(Enum):
    ONLINE = 1
    """
    Services that are always online. If these fail, we fail the workflow and
    log this.
    """
    
    INTERMITTENT_ERROR_CODE = 2
    """
    Services that error often (not go down!)
    these will be logged when they fail, but we continue with the cached option.
    """

    # (we choose to do this over arbitrary lambdas because its nicer. For now.)
    INTERMITTENT_HTML = 3
    """
    Like INTERMITTENT_ERROR_CODE, but errors when HTML is returned.
    """

@dataclass
class CacheItem:
    """
    Class for differentriating between offline and online items in a cache.

    NOTE: If cached is "", we assume that online is a Google Drive URL (for cases where there is no
    remaining online data source.)
    """

    name: str
    """The display name of the artifact, used for human-printing."""
    cached: str
    online: Optional[Service] = None
    status: OnlineStatus = OnlineStatus.ONLINE
    """How much to care about errors from downloading the online file."""

    @classmethod
    @warnings.deprecated("Pending for removal after the CONTRIBUTING guide is updated.")
    def cache_only(cls, name: str, cached: str) -> "CacheItem":
        """Wrapper method to explicitly declare a CacheItem as cached only."""
        return cls(name=name, cached=cached, online=None)

    def download(self, output: str | PathLike):
        print(f"Fetching {self.name}...")

        with NamedTemporaryFile() as cached_file:
            print(f"Downloading cache {self.cached}...")
            gdown.download(self.cached, cached_file)

            if self.online is None: return

            print(f"Downloading {self.online}...")
            with self.online.download(output) as response:

                print("Checking that downloaded artifact matches with cached artifact...")
                if filecmp.cmp(output, cached_file.name): return # It does!

                # For debug purposes, we allow the output artifact to be viewed in some kind of temporary folder.
                debug_file_path = Path(NamedTemporaryFile(prefix="spras-benchmarking-debug-artifact", delete=False).name)
                # (and we pedantically use this over Path#rename since temporary directories can be mounted to a different file system.)
                shutil.move(output, debug_file_path)
                if (self.status == OnlineStatus.INTERMITTENT_ERROR_CODE and not response.ok) \
                    or (self.status == OnlineStatus.INTERMITTENT_HTML and Path(debug_file_path).read_text().strip().startswith("<!DOCTYPE html>")):
                    warnings.warn(f"Online url {self.online} erroring with status code {response.status_code}. See {debug_file_path} for the online output. Using the cached file instead...")
                    # Back up to the cached_file
                    shutil.move(cached_file.name, output)
                else:
                    raise RuntimeError(f"Cached and online files did not match with status code {response.status_code}! See {debug_file_path} for the online output.")


CacheDirectory = dict[str, Union[CacheItem, "CacheDirectory"]]

# An *unversioned* directory list.
directory: CacheDirectory = {
    "STRING": {
        "9606": {
            "9606.protein.links.full.txt.gz": CacheItem(
                name="STRING 9606 full links",
                cached="https://drive.google.com/uc?id=13tE_-A6g7McZs_lZGz9As7iE-5cBFvqE",
                online=Service("http://stringdb-downloads.org/download/protein.links.full.v12.0/9606.protein.links.full.v12.0.txt.gz"),
            ),
            "9606.protein.aliases.txt.gz": CacheItem(
                name="STRING 9606 protein aliases",
                cached="https://drive.google.com/uc?id=1IWrQeTVCcw1A-jDk-4YiReWLnwP0S9bY",
                online=Service("https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"),
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
                online=Service("https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names&format=tsv&query=%28*%29+AND+%28reviewed%3Atrue%29+AND+%28model_organism%3A9606%29"),
            ),
            # idmapping FTP files. See the associated README:
            # https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
            "HUMAN_9606_idmapping_selected.tab.gz": CacheItem(
                name="UniProt 9606 ID external database mapping",
                cached="https://drive.google.com/uc?id=1Oysa5COq31H771rVeyrs-6KFhE3VJqoX",
                online=Service("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz"),
            ),
            "HUMAN_9606_idmapping.dat.gz": CacheItem(
                name="UniProt 9606 internal id mapping",
                cached="https://drive.google.com/uc?id=1lGxrx_kGyNdupwIOUXzfIZScc7rQKP-O",
                online=Service("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz"),
            ),
        }
    },
    "DISEASES": {
        # Instead of going through https://unmtid-shinyapps.net/shiny/tiga/, we use their
        # archived files directory instead.
        "tiga_gene-trait_stats.tsv": CacheItem(
            name="TIGA data",
            cached="https://drive.google.com/uc?id=114qyuNDy4qdmYDHHJAW-yBeTxcGTDUnK",
            online=Service("https://unmtid-dbs.net/download/TIGA/20250916/tiga_gene-trait_stats.tsv"),
        ),
        "HumanDO.tsv": CacheItem(
            name="Disease ontology data",
            cached="https://drive.google.com/uc?id=1lfB1DGJgrXTxP_50L6gGu_Nq6OyDjiIi",
            online=Service("https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/016a4ec33d1a1508d669650086cd92ccebe138e6/DOreports/HumanDO.tsv"),
        ),
        "human_disease_textmining_filtered.tsv": CacheItem(
            name="DISEASES textmining channel",
            cached="https://drive.google.com/uc?id=1vD8KbT9sk04VEJx9r3_LglCTGYJdhN0D",
            online=Service("https://download.jensenlab.org/human_disease_textmining_filtered.tsv"),
        ),
        "human_disease_knowledge_filtered.tsv": CacheItem(
            name="DISEASES knowledge channel",
            cached="https://drive.google.com/uc?id=1qGUnjVwF9-8p5xvp8_6CfVsbMSM_wkld",
            online=Service("https://download.jensenlab.org/human_disease_knowledge_filtered.tsv"),
        ),
    },
    "BioMart": {
        "ensg-ensp.tsv": CacheItem(
            name="BioMart ENSG <-> ENSP mapping",
            cached="https://drive.google.com/uc?id=1-gPrDoluXIGydzWKjWEnW-nWhYu3YkHL",
            online=fetch_biomart_service((dir_path / "biomart" / "ensg-ensp.xml").read_text()),
        )
    },
    "DepMap": {
        "OmicsProfiles.csv": CacheItem(
            name="DepMap omics metadata",
            cached="https://drive.google.com/uc?id=1i54aKfO0Ci2QKLTNJnuQ_jgGhH4c9rTL",
            online=Service("https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2F2025-05-01-master-mapping-table-28c2.12%2Fpublic_release_date.2025-05-01.master_mapping_table.csv&dl_name=OmicsProfiles.csv&bucket=depmap-external-downloads"),
        ),
        "CRISPRGeneDependency.csv": CacheItem(
            name="DepMap gene dependency probability estimates",
            cached="https://drive.google.com/uc?id=122rWNqT_u3M7B_11WYZMtOLiPbBykkaz",
            online=Service("https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2F25q2-public-557c.3%2FCRISPRGeneDependency.csv&dl_name=CRISPRGeneDependency.csv&bucket=depmap-external-downloads"),
        ),
        "OmicsSomaticMutationsMatrixDamaging.csv": CacheItem(
            name="DepMap genotyped matrix",
            cached="https://drive.google.com/uc?id=1W7N2H0Qi7NwmTmNChcwa2ZZ4WxAuz-Xh",
            online=Service("https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2Fpublic-25q2-c5ef.87%2FOmicsSomaticMutationsMatrixDamaging.csv&dl_name=OmicsSomaticMutationsMatrixDamaging.csv&bucket=depmap-external-downloads"),
        ),
        "OmicsExpressionProteinCodingGenesTPMLogp1.csv": CacheItem(
            name="DepMap model-level TPMs",
            cached="https://drive.google.com/uc?id=1P0m88eXJ8GPdru8h9oOcHPeXKU7ljIrP",
            online=Service("https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2Fpublic-25q2-c5ef.73%2FOmicsExpressionProteinCodingGenesTPMLogp1.csv&dl_name=OmicsExpressionProteinCodingGenesTPMLogp1.csv&bucket=depmap-external-downloads"),
        ),
        "OmicsCNGeneWGS.csv": CacheItem(
            name="DepMap gene-level copy number data",
            cached="https://drive.google.com/uc?id=1TPp3cfK7OZUrftucr3fLO-krXSQAA6Ub",
            online=Service("https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2Fpublic-25q2-c5ef.104%2FOmicsCNGeneWGS.csv&dl_name=OmicsCNGeneWGS.csv&bucket=depmap-external-downloads"),
        ),
    },
    "iRefIndex": {
        # This can also be obtained from the SPRAS repo
        # (https://github.com/Reed-CompBio/spras/blob/b5d7a2499afa8eab14c60ce0f99fa7e8a23a2c64/input/phosphosite-irefindex13.0-uniprot.txt).
        # iRefIndex has been down for quite some time, so this is only from the cache.
        "phosphosite-irefindex13.0-uniprot.txt": CacheItem.cache_only(
            name="iRefIndex v13.0 UniProt interactome",
            cached="https://drive.google.com/uc?id=1fQ8Z3FjEwUseEtsExO723zj7mAAtdomo"
        )
    },
    "OsmoticStress": {
        "yeast_pcsf_network.sif": CacheItem.cache_only(
            # In the paper https://doi.org/10.1016/j.celrep.2018.08.085
            name="Case Study Edge Results, from Supplementary Data 3",
            cached="https://drive.google.com/uc?id=1Agte0Aezext-8jLhGP4GmaF3tS7gHX-h"
        ),
        # The following files are from https://github.com/gitter-lab/osmotic-stress
        "prizes.txt": CacheItem(
            name="Osmotic Stress Prizes",
            online=Service("https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/Input%20Data/prizes.txt"),
            cached="https://drive.google.com/uc?id=16WDQs0Vjv6rI12-hbifsbnpH31jMGhJg"
        ),
        "ChasmanNetwork-DirUndir.txt": CacheItem(
            name="Network Input",
            online=Service("https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/Input%20Data/ChasmanNetwork-DirUndir.txt"),
            cached="https://drive.google.com/uc?id=1qYXPaWcPU72YYME7NaBzD7thYCHRzrLH"
        ),
        "dummy.txt": CacheItem(
            name="Dummy Nodes File",
            online=Service("https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/Input%20Data/dummy.txt"),
            cached="https://drive.google.com/uc?id=1dsFIhBrIEahggg0JPxw64JwS51pKxoQU"
        ),
        "_edgeFreq.eda ": CacheItem(
            name="Case Study Omics Integrator Edge Frequencies",
            online=Service("https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/Notebooks/Forest-TPS/_edgeFreq.eda"),
            cached="https://drive.google.com/uc?id=1M_rxEzUCo_EVuFyM47OEH2J-4LB3eeCR"
        ),
        "goldStandardUnionDetailed.txt": CacheItem(
            name="Gold Standard Reference Pathways",
            online=Service("https://raw.githubusercontent.com/gitter-lab/osmotic-stress/refs/heads/master/data/evaluation/goldStandardUnionDetailed.txt"),
            cached="https://drive.google.com/uc?id=1-_zF9oKFCNmJbDCC2vq8OM17HJw80s2T"
        ),
    },
    "Surfaceome": {
        "table_S3_surfaceome.xlsx": CacheItem(
            name="Human surfaceome",
            online=Service("http://wlab.ethz.ch/surfaceome/table_S3_surfaceome.xlsx"),
            cached="https://docs.google.com/uc?id=1cBXYbDnAJVet0lv3BRrizV5FuqfMbBr0"
        )
    },
    "TranscriptionFactors": {
        "Homo_sapiens_TF.tsv": CacheItem.cache_only(
            name="Human transcription factors",
            # This server has anti-bot protection, so to respect their wishes, we don't download from the server.
            # The original URL is https://guolab.wchscu.cn/AnimalTFDB4_static/download/TF_list_final/Homo_sapiens_TF,
            # which is accessible from https://guolab.wchscu.cn/AnimalTFDB4//#/Download -> Homo sapiens
            # (also under the Internet Archive as of Feb 2nd, 2026. If the original artifact disappears, the drive link below should suffice.)
            cached="https://drive.google.com/uc?id=1fVi18GpudUlquRPHgUJl3H1jy54gO-uz",
        )
    },
    "PathwayCommons": {
        # TODO: all of these share the same common URL: can we make this API a little nicer?
        "PANTHER": {
            "Apoptosis_signaling_pathway.txt": CacheItem(
                name="Apoptosis Signaling Pathway",
                cached="https://drive.google.com/uc?id=1BPcnvqHrGMQeX4oQx2ow3OribgPxzwhG",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00006"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "B_cell_activation.txt": CacheItem(
                name="B cell activation",
                cached="https://drive.google.com/uc?id=1iWcb5AfdobGncRB6xQ6T5qunXzb6Gxd-",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00010"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "Beta3_adrenergic_receptor_signaling_pathway.txt": CacheItem(
                name="Beta3_adrenergic_receptor_signaling_pathway",
                cached="https://drive.google.com/uc?id=1jrJzrDvhDAs818wYjQ_dm1irOz8Bv4lk",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP04379"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "Cadherin_signaling_pathway.txt": CacheItem(
                name="Cadherin signaling pathway",
                cached="https://drive.google.com/uc?id=14Of-6mwIpul_QciyJ-Xb9f7t-IrVcIna",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00012"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "Fas_signaling_pathway.txt": CacheItem(
                name="FAS signaling_pathway",
                cached="https://drive.google.com/uc?id=121cHJf0ZtglQHvy9xuEpYSjwBbJV9Fju",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00020"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "FGF_signaling_pathway.txt": CacheItem(
                name="FGF signaling pathway",
                cached="https://drive.google.com/uc?id=1PIiWK1-ImXE1YHdDh1hGUVB01Ye8brQg",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00021"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "Hedgehog_signaling_pathway.txt": CacheItem(
                name="Hedgehog signaling pathway",
                cached="https://drive.google.com/uc?id=1i7HKn4nlJQcaXUDXpbpDFBxbkBXZC0xQ",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00025"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "Insulin_IGF_pathway_protein_kinase_B_signaling_cascade.txt": CacheItem(
                name="Insulin/IGF pathway-protein kinase B signaling cascade",
                cached="https://drive.google.com/uc?id=1Xkxcm0ngrE8otau9ccyPeCg7KZUdhJf7",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00033"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "Interferon_gamma_signaling_pathway.txt": CacheItem(
                name="Interferon-gamma signaling pathway",
                cached="https://drive.google.com/uc?id=1aPqi0A5ZIOA5kKELVUI_NvC8taiHll5z",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00035"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "Interleukin_signaling_pathway.txt": CacheItem(
                name="Interleukin signaling pathway",
                cached="https://drive.google.com/uc?id=1IOv14pRJ8aN9LRnkZ4BQXf3QGUAashku",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00036"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "JAK_STAT_signaling_pathway.txt": CacheItem(
                name="JAK/STAT signaling pathway",
                cached="https://drive.google.com/uc?id=1QzMEMUZzeoxUYZZRGcm6Al_HzH6pmwED",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00038"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "Nicotinic_acetylcholine_receptor_signaling_pathway.txt": CacheItem(
                name="Nicotinic acetylcholine receptor signaling pathway",
                cached="https://drive.google.com/uc?id=1SdnKr4TthfmZWgMA_FOlTmf-EEpNsdzx",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00044"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "Notch_signaling_pathway.txt": CacheItem(
                name="Notch signaling pathway",
                cached="https://drive.google.com/uc?id=1qfyxuc1EomOKGRyI7QyQ7LUUhLPZytz5",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00045"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "PDGF_signaling_pathway.txt": CacheItem(
                name="PDGF signaling pathway",
                cached="https://drive.google.com/uc?id=1A9hl340XKnZeNfd3hiiX7lxOVV94lQ5s",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00047"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "Ras_pathway.txt": CacheItem(
                name="Ras pathway",
                cached="https://drive.google.com/uc?id=1wNizL5wDh48E5YxHcZjURa9UeKMONrgr",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP04393"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "T_cell_activation.txt": CacheItem(
                name="T cell activation",
                cached="https://drive.google.com/uc?id=1t5G_jN8QSOiVceQGAmKvbYebkV1G5oJy",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00053"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "Toll_receptor_signaling_pathway.txt": CacheItem(
                name="Toll receptor signaling pathway",
                cached="https://drive.google.com/uc?id=1nFix8mMvuU_Vu9tExwgaS279nynqM_Oo",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00054"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "VEGF_signaling_pathway.txt": CacheItem(
                name="VEGF signaling pathway",
                cached="https://drive.google.com/uc?id=1W1G0TmA6-JLF9pIZD0TR4w95IwG2IALs",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00056"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
            "Wnt_signaling_pathway.txt": CacheItem(
                name="Wnt signaling pathway",
                cached="https://drive.google.com/uc?id=1diaacbik5hcA9Fo7vMXFAP_wXRe0xCLB",
                online=Service("https://www.pathwaycommons.org/pc2/get?format=TXT&uri=https%3A%2F%2Fidentifiers.org%2Fpanther.pathway%3AP00057"),
                status=OnlineStatus.INTERMITTENT_HTML
            ),
        }
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
