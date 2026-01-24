# cache

Handles artifact fetching and cache. This folder has:

- `Snakefile` which only contains a function used for producing fetching rules.
- `directory.py`, the actual location of file URLs and their cached counterparts.
- `cli.py`, a utility for manually fetching specific URLs from `directory.py`.
- `util.py`, an internal file for use by the other files above.
