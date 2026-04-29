"""
This is how spras-benchmarking handles artifact downloading and caching.
`cache.link` should be used specifically inside `Snakefile`
"""

from dataclasses import dataclass
from pathlib import Path
import pickle
from urllib.parse import quote_plus
from typing import Optional, Union

from cache.util import PostProcessAction
from cache.directory import CacheItem, get_cache_item


__all__ = ["FetchConfig", "link"]

dir_path = Path(__file__).parent.resolve()
artifacts_dir = dir_path / "artifacts"


@dataclass(frozen=True)
class FetchConfig:
    """
    What `directive` (a CacheItem or a location in `./directory.py`)
    should be fetched, and how it should be fetched.
    """

    directive: Union[CacheItem, tuple[str, ...]]
    post_process: Optional[PostProcessAction] = None
    # NOTE: uncompress only unzips `.gz` files. TODO: add support for .zip files.


def stringify_tuple_directive(directive: tuple[str, ...]) -> str:
    return quote_plus("/".join(directive))


def stringify_config(directive: Union[CacheItem, tuple[str, ...]]) -> str:
    """
    Directives aren't friendly for Snakemake rules: we use urllib to make a URL-safe string, as Snakemake accepts
    such strings as rule names. They aren't the most readable, but they suffice for logging purposes.
    """
    return quote_plus(directive.name if isinstance(directive, CacheItem) else "/".join(directive))


def add_suffix(path: Path, suffix: str):
    """`file.suffix` -> `file.suffix.suffix2`."""
    return path.with_suffix(path.suffix + suffix)


def metadata_has_expired(cache_item: CacheItem, output: Path) -> bool:
    """
    Check if the artifact metadata (the serialized `cache_item`) associated with a `cache_item` has expired (i.e. changed or deleted).
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

    # metadata does not match the `cache_item`: re-retrieve the item
    if old_cache_item != cache_item:
        with open(metadata_file, "wb") as f:
            pickle.dump(cache_item, f)
        return True

    # metadata hasn't changed and already existed: this hasn't expired
    return False


def link_with_cache_item(output: Path, cache_item: CacheItem, post_process: Optional[PostProcessAction] = None):
    """
    Intermediary function for `link`.
    This does almost all of what `link` is characterized to do in its documentation,
    except for doing symlinking.
    """
    # If `uncompress` is `True`, we make
    # `output` our 'compressed output.'
    uncompressed_output = output
    if post_process is not None:
        output = add_suffix(output, ".unprocesed")

    # Re-download if the file doesn't exist or the directive has expired.
    # Note that we check for expiration first to trigger metadata creation.
    if metadata_has_expired(cache_item, output) or not output.exists():
        output.unlink(missing_ok=True)
        cache_item.download(output)

    if post_process is not None:
        processed_artifact_path = add_suffix(output, ".processed")
        processed_artifact_path.unlink(missing_ok=True)
        post_process.run_action(output, uncompressed_output)


def link(output: str, config: FetchConfig):
    """
    Links output files from the provided `config` to the provided `output`.

    @param output: The output file to write to.
    @param config: The config to fetch from.

    For example,

    ```python
    link("output/ensg-ensp.tsv", FetchConfig(("BioMart", "ensg-ensp.tsv")))
    ```

    Will download and check BioMart's cache for ENSG-ENSP mapping, then:
    - If `config.directive` is a `CacheItem`, we write the file directly to `output`.
    - Otherwise, we symlink the cached output (lying in `artifacts_dir`) with the desired `output`
    to avoid unnecessary file duplication.

    In the above example's case, since the configuration points to an entry in `directory.py`,
    we would do the latter symlinking option, avoiding querying the BioMart API in case this data
    is requested again.

    This function wraps around `link_with_cache_item` and handles symlinking
    depending on the type of config.directive.
    """
    # TODO: there is most likely a nicer way to design this.

    if isinstance(config.directive, CacheItem):
        link_with_cache_item(Path(output), config.directive, config.post_process)
    else:
        artifacts_dir.mkdir(exist_ok=True)
        artifact_name = stringify_tuple_directive(config.directive)
        artifact_output = artifacts_dir / artifact_name

        link_with_cache_item(artifact_output, get_cache_item(config.directive), config.post_process)

        Path(output).symlink_to(artifact_output)
