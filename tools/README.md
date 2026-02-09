# Dataset Processing Tools

This includes common tools for doing dataset processing, which take in SPRAS-compatible file formats. This currently includes:

- `trim.py`: this may be included in SPRAS later, but this contains utilities for trimming a gold standard with its respective interactome,
and the gold standard data with the interactome and the gold standard itself.
- `sample.py`: this samples an interactome and downstream samples the gold standard, preserving a percentage of the associated data in the largest
connected component of the gold standard. _These tools require a gold standard_
