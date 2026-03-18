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
    # Either-or
    pinned=Service("..."),
    unpinned=Service("..."),
),
```

When a file is requested, `cached`, `pinned`, and `unpinned` are all downloaded:
- If the URLs linking to `pinned` and `unpinned` do not succeed (i.e. do not return a 2XX status code), we fail.
- If the URL linking to `pinned` does not match `cached`, we fail.
- If the URL linking to `unpinned` does not match `cached`, we warn that the data needs updating.

## Layout

This folder has:
- `Snakefile` which only contains a function used for producing fetching rules.
- `directory.py`, the actual location of file URLs and their cached counterparts.
- `cli.py`, a utility for manually fetching specific URLs from `directory.py`.
- `util.py`, an internal file for use by the other files above.
