"""
This is how spras-benchmarking handles artifact caching. `cache` should be used specifically inside `Snakefile`
"""

from cache.util import uncompress as uncompress_file
from cache.directory import get_cache_item
from pathlib import Path
import os
from urllib.parse import quote_plus
import pickle

__all__ = ["link"]

dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
artifacts_dir = dir_path / "artifacts"

def get_artifact_name(directive: list[str]) -> str:
    return quote_plus("/".join(directive))

def has_expired(directive: list[str]) -> bool:
    """
    Check if the artifact metadata associated with a directive has expired.
    Avoids re-downloading the artifact if nothing has changed.
    """
    artifact_name = get_artifact_name(directive)
    cache_item = get_cache_item(directive)

    metadata_dir = artifacts_dir / 'metadata'
    metadata_dir.mkdir(exist_ok=True)
    metadata_file = (artifacts_dir / 'metadata' / artifact_name).with_suffix((artifacts_dir / artifact_name).suffix + '.metadata')

    # metadata never existed: we need to retrieve the new file
    if not metadata_file.exists():
        with open(metadata_file, 'wb') as f:
            pickle.dump(cache_item, f)
        return True

    old_cache_item = None
    with open(metadata_file, 'rb') as f:
        old_cache_item = pickle.load(f)

    # metadata expired: re-retrieve the item
    if old_cache_item != cache_item:
        with open(metadata_file, 'wb') as f:
            pickle.dump(cache_item, f)
        return True

    # metadata hasn't changed and already existed: this hasn't expired
    return False

def link(output: str, directive: list[str], uncompress=False):
    """
    Links output files from cache.directory directives.
    For example,

    ```py
    link("output/ensg-ensp.tsv", ["BioMart", "ensg-ensp.tsv"])
    ```

    would download and check BioMart's cache for ENSG-ENSP mapping, then symlink the cached output
    (lying somewhere in the cache folder) with the desired `output`.
    """

    artifacts_dir.mkdir(exist_ok=True)

    artifact_name = get_artifact_name(directive)

    Path(output).unlink(missing_ok=True)

    # Re-download if the directive has expired.
    cache_item = get_cache_item(directive)
    if has_expired(directive):
        (artifacts_dir / artifact_name).unlink(missing_ok=True)
        cache_item.download(artifacts_dir / artifact_name)

    if uncompress:
        uncompressed_artifact_path = Path(str(artifacts_dir / artifact_name) + '.uncompressed')
        uncompressed_artifact_path.unlink(missing_ok=True)
        uncompress_file(artifacts_dir / artifact_name, uncompressed_artifact_path)
        Path(output).symlink_to(uncompressed_artifact_path)
    else:
        Path(output).symlink_to(artifacts_dir / artifact_name)
