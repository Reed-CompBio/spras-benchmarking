"""
This is how spras-benchmarking handles artifact caching. `cache` should be used specifically inside `Snakefile`
"""

from dataclasses import dataclass
from typing import Union
from cache.util import uncompress as uncompress_file
from cache.directory import CacheItem, get_cache_item
from pathlib import Path
import os
from urllib.parse import quote_plus
import pickle

__all__ = ["FetchConfig", "link"]

dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
artifacts_dir = dir_path / "artifacts"

@dataclass(frozen=True)
class FetchConfig:
    directive: Union[CacheItem, tuple[str, ...]]
    uncompress: bool = False

def get_artifact_name(directive: tuple[str, ...]) -> str:
    return quote_plus("/".join(directive))

def add_suffix(path: Path, suffix: str):
    return path.with_suffix(path.suffix + suffix)

def has_expired(
        cache_item: CacheItem,
        output: Path
) -> bool:
    """
    Check if the artifact metadata associated with a directive has expired.
    Avoids re-downloading the artifact if nothing has changed.
    """

    metadata_file = add_suffix(output, ".metadata")

    # metadata never existed: we need to retrieve the new file
    if not metadata_file.exists():
        with open(metadata_file, "wb") as f:
            pickle.dump(cache_item, f)
        return True

    old_cache_item = None
    with open(metadata_file, "rb") as f:
        old_cache_item = pickle.load(f)

    # metadata expired: re-retrieve the item
    if old_cache_item != cache_item:
        with open(metadata_file, "wb") as f:
            pickle.dump(cache_item, f)
        return True

    # metadata hasn't changed and already existed: this hasn't expired
    return False

def link_with_cache_item(
    output: Path,
    cache_item: CacheItem,
    uncompress: bool = False
):
    """
    Intermediary function for `link`.
    This does almost all of what `link` is characterized to do in its documentation,
    except for doing symlinking.
    """
    # If `uncompress` is `True`, we make
    # `output` our 'compressed output.'
    uncompressed_output = output
    if uncompress:
        output = add_suffix(output, ".compressed")

    # Re-download if the file doesn't exist or the directive has expired.
    # Note that we check for expiration first to trigger metadata creation.
    if has_expired(cache_item, output) or not output.exists():
        output.unlink(missing_ok=True)
        cache_item.download(output)

    if uncompress:
        uncompressed_artifact_path = add_suffix(output, ".uncompressed")
        uncompressed_artifact_path.unlink(missing_ok=True)
        uncompress_file(output, uncompressed_output)

def link(
        output: str,
        config: FetchConfig
):
    """
    Links output files from cache.directory directives.
    For example,

    ```py
    link("output/ensg-ensp.tsv", FetchConfig(("BioMart", "ensg-ensp.tsv")))
    ```

    would download and check BioMart's cache for ENSG-ENSP mapping, then:
    - If `config.directive` is a `CacheItem`, we write the file directly to `output`.
    - Otherwise, we symlink the cached output (lying somewhere in the cache folder) with the desired `output`
    to avoid file duplication.

    This function wraps around link_with_cache_item and handles symlinking
    depending on the type of config.directive.
    TODO: most likely a nicer way to design this.
    """

    if isinstance(config.directive, CacheItem):
        link_with_cache_item(
            Path(output),
            config.directive,
            config.uncompress
        )
    else:
        artifacts_dir.mkdir(exist_ok=True)
        artifact_name = get_artifact_name(config.directive)
        artifact_output = artifacts_dir / artifact_name

        link_with_cache_item(
            artifact_output,
            get_cache_item(config.directive),
            config.uncompress
        )

        Path(output).symlink_to(artifact_output)

