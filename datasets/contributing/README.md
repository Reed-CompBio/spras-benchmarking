# Contributing Guide dataset

**This is an artificial dataset** for how to make datasets.

This comes with a `raw_generation.py` script, which produces the associated raw data, where the gold standard is `k` paths of length `n` with
Erdős-Rényi edges, such that the sources and targets come from the start and ends of each path. The background interactome is the gold standard with
more edge and node noise. This is not a topologically-accurate emulation of (signaling) pathways, but it suffices to trick most pathway reconstruction
algorithms.

This does not cover the (very common!) task of ID mapping, as this can vary constantly between datasets.
