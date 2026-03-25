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

When a file is requested, `cached`, `pinned`, and `unpinned` are all downloaded. `cached` is the link to the underlying file that we store,
`pinned` is the link to an arbitrary online service containing a versioned file that never changes (we use this to check for uptime),
and `unpinned` is the link to the an arbitrary online service containing an unversioned file we use to check for updates.

We characterize them as follows:
- If the URLs linking to `pinned` and `unpinned` do not succeed (i.e. do not return a 2XX status code), we fail.
- If the URL linking to `pinned` does not match `cached`, we fail.
- If the URL linking to `unpinned` does not match `cached`, we warn that the data needs updating. The data itself will not automatically update.

Specifically, `unpinned` links to file URLs that constantly update, `pinned` does otherwise, and `cached` links to our
own copy of the data that should match with the `unpinned` and `pinned` URLs. We prefer to have both `pinned` and `unpinned` URLs, but
there are many situations where the `pinned` URL is not available (e.g. the queried service has no versioning), or the `unpinned` URL is not available
(e.g. the queried service only has versioning).

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
    "raw/9606.protein.links.full.txt": FetchConfig(("STRING", "9606", "9606.protein.links.full.txt.gz"), uncompress=True),
})
```

would produce a Snakemake rule whose output is `raw/9606.protein.links.full.txt`, and would look under `directory.py` by traversing
the `directory` dictionary, going to `STRING` then `9606` then `9606.protein.links.full.txt.gz`, where the `FetchConfig`
asks the inner rule to uncompress the file before saving it under `raw/9606.protein.links.full.txt`.

Semantically, this is equivalent to:

```py
produce_fetch_rules({
    "raw/9606.protein.links.full.txt": FetchConfig(CacheItem(
        name="STRING 9606 full protein links",
        cached="...",
        pinned="...",
    ), uncompress=True),
})
```

However, the former option, since it uses items in `directory.py`, saves the file to a cached directory under `cache/artifacts`.
The latter saves the file to a dataset-specific folder for dataset `Snakefile`s: that is, if you have a file
that's used across multiple datasets, add it to `directory.py`!

## Implementation details

### `.metadata`

All cached files come with an associated `.metadata`: usually, this would be controlled with Snakemake, but since this system lives
outside of the purview of `Snakemake`, we instead track file data with an associated `.metadata` file, which preserves information
about where the file came from, and when it was created, to re-fetch files if any of that associated data changes.
This is controlled under `__init__.py`.

### Loggers

Later on, our use of `loguru` will be logged to let maintainers know when data sources are outdated.

### Layout

This folder has:
- `Snakefile` which only contains a function used for producing fetching rules.
- `directory.py`, where named `CacheItem`s are stored, as well as the code that defines
the schema (including `CacheItem`) for the rest of this directory.
- `cli.py`, a debugging utility for manually fetching specific data from `directory.py`.
- `util.py`, an internal file for use by the other files above.
- `__init__.py`, which acts as an intermediary between `Snakefile` and `directory.py`, providing utilities for handling file metadata.
