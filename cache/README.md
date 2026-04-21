# Cache

Handles artifact fetching and cache. The point of this is to [for the duration that SPRAS is maintained] prevent any kind of
data rot, to ensure that continuous benchmarking is encouraged to use the latest available data.

During benchmarking runs, data is fetched from all provided URLs in `directory.py`, where we get the most current version of data,
and compare it to our cached data to check if the data has changed at all.

All entries are provided with this template:

```py
"file-name.ext": CacheItem(
    name="Short File Description",
    cached="https://drive.google.com/uc?id=...",
    # These two are optional, but heavily encouraged to be included.
    pinned="...",
    unpinned="...",
),
```

When a file is requested, `cached`, `pinned`, and `unpinned` are all downloaded. `cached` is the link to the underlying file that we store in Drive,
`pinned` is the link to an arbitrary online service containing a versioned static file (we use this to check for uptime and dually to document
where a file was originally retrieved), and `unpinned` is the link to the an arbitrary online service containing an unversioned
(which has the capability of becoming newer) file we use to check for updates.

We characterize them as follows:
- If the URL linking to `pinned` does not match `cached`, we fail, as this means that our versioned file updated/errored in some way.
- If the URL linking to `unpinned` does not match `cached`, we warn that the data needs updating. The data itself will not automatically update.

Specifically, `unpinned` links to file URLs that constantly update, `pinned` does otherwise, and `cached` links to our
own copy of the data that must match with  the `pinned` URL, and is encouraged to match with the `unpinned` URL to prevent our own copy of the data
from drifting with respect to its version. We prefer to have both `pinned` and `unpinned` URLs, but
there are many situations where the `pinned` URL is not available (e.g. the queried service has no versioning), or the `unpinned` URL is not available
(e.g. the queried service only has versioning).

If `cache` doesn't match `pinned`, this usually indicates that the service is down: we don't have a way to handle temporary outages at the moment,
but permanent outages should remove references to `pinned` entirely, noting that the linked service is down forever for some reason.

## Google Drive

We currently use Google Drive to store raw data. The hope is to move to [OSDF](https://osg-htc.org/services/osdf).
We have been running into the occasional ratelimiting issue, which may become more of a problem in the future.

## Snakemake

We also provide a `Snakefile`, which contains dataset fetching functions that can be imported in dataset-specific Snakefiles through:

```py
include: "../../cache/Snakefile"
```

This imports a function `produce_fetch_rules`, which takes in a dictionary where the keys are file names,
and the values are either entries in `directory.py`, or `CacheItem`s themselves. For example,

```py
produce_fetch_rules({
    "raw/9606.protein.links.full.txt": FetchConfig(
        ("STRING", "v12", "9606", "9606.protein.links.full.txt.gz"),
        post_process=PostProcessAction.UNCOMPRESS_GZ
    ),
})
```

would produce a Snakemake rule whose output is `raw/9606.protein.links.full.txt`, and would look under `directory.py` by traversing
the `directory` dictionary, going to `STRING` then `v12` then `9606` then `9606.protein.links.full.txt.gz`, where the `FetchConfig`
asks the inner rule to uncompress the .gz file before saving it under `raw/9606.protein.links.full.txt`.

> [!NOTE]
> `FetchConfig` is optional, and is only used to specify the `post_process` option.

Semantically, this is equivalent to:

```py
produce_fetch_rules({
    "raw/9606.protein.links.full.txt": FetchConfig(CacheItem(
        name="STRING 9606 full protein links",
        cached="...",
        pinned="...",
    ), post_process=PostProcessAction.UNCOMPRESS_GZ),
})
```

However, the former option, since it uses items in `directory.py`, saves the file to a cached directory under `cache/artifacts`.
The latter saves the file to a dataset-specific folder for dataset `Snakefile`s: that is, if you have a file
that's used across multiple datasets, add it to `directory.py`. Otherwise, if you have a file specifically used for a dataset, keep it under that
respective dataset.

Some datasets also provide an `unversioned` tag, which would point to the latest version under `directory`. For example, at the time of writing
(when `unversioned` pointed to `v12`), the following are equivalent:

```py
FetchConfig(("STRING", "v12", "9606", "9606.protein.links.full.txt.gz"), post_process=PostProcessAction.UNCOMPRESS_GZ)
FetchConfig(("STRING", "unversioned", "9606", "9606.protein.links.full.txt.gz"), post_process=PostProcessAction.UNCOMPRESS_GZ)
```

## Implementation details

### `.metadata`

All cached files come with an associated `.metadata`: usually, this would be controlled with Snakemake, but since this system lives
outside of the purview of `Snakemake`, we instead track file data with an associated `.metadata` file, which preserves information
about where the file came from, and when it was created, to re-fetch files if any of that associated data changes.
This is controlled under `link.py`.

### Loggers

Later on, our use of `loguru` will be logged to let maintainers know when data sources are outdated.

### Layout

This folder has:
- `Snakefile` which only contains a function used for producing fetching rules.
- `directory.py`, where named `CacheItem`s are stored, as well as the code that defines
the schema (including `CacheItem`) for the rest of this directory.
- `cli.py`, a debugging utility for manually fetching specific data from `directory.py`.
- `util.py`, an internal file for use by the other files above.
- `link.py`, which acts as an intermediary between `Snakefile` and `directory.py`, providing utilities for handling file metadata.
