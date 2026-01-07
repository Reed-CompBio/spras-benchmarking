"""
This is how spras-benchmarking handles artifact caching. `cache` should be used specifically inside `Snakefile`
"""

from cache.directory import get_cache_item
from pathlib import Path
import os
from urllib.parse import quote_plus

__all__ = ["link"]

dir_path = Path(os.path.dirname(os.path.realpath(__file__)))

def link(output: str, directive: list[str]):
    """
    Links output files from cache.directory directives.
    For example,
    
    ```py
    link("output/ensg-ensp.tsv", ["BioMart", "ensg-ensp.tsv"])
    ```

    would download and check BioMart's cache for ENSG-ENSP mapping, then symlink the cached output
    (lying somewhere in the cache folder) with the desired `output`.
    """

    artifacts_dir = dir_path / "artifacts"
    artifacts_dir.mkdir(exist_ok=True)

    artifact_name = quote_plus("/".join(directive))

    get_cache_item(directive).download(artifacts_dir / artifact_name)
    (artifacts_dir / artifact_name).symlink_to(output)
