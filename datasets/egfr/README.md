# EGFR

EGFR dataset. This dataset does a lot less processing for raw files, and is mainly focused on creating the new STRING-based interactome.

## Overview

This produces two sets of files: one based on the iRefIndex/PhosphoSite directed interactome of closed-source origin based off of UniProt identifiers, and another one based off of the more updated though undirected STRING interactome.

## Data

_See `Snakefile` for specific file locations._

The score data (`egfr-prizes.txt`), gold standard nodes `eight-egfr-reference-all.txt`, and the (now-deprecated) manually edited `iRefIndex`-based interactome are all from [_Synthesizing Signaling Pathways from Temporal Phosphoproteomic Data_](https://doi.org/10.1016/j.celrep.2018.08.085).

We also use the StringDB human interactome and UniProt mapping files. See `cache/directory.py` for more info on these.

## Scripts

- `process_prizes.py`: produces a `prizes-uniprot.txt` from
[egfr-prizes.txt](https://raw.githubusercontent.com/gitter-lab/tps/refs/heads/master/data/pcsf/egfr-prizes.txt),
trimming psuedonodes and manually injecting the `EGF_HUMAN` receptor as a dummy node for OmicsIntegrator1.
- `process_interactome.py`: Produces the STRING `interactome.tsv` file from the STRING links file. Note that the `phosphosite-irefindex13.0-uniprot.txt` is a magic (as in with closed-source origin) directed interactome produced with a combination of the now archived iRefIndex v13 interactome with extra PhosphoSite-provided nodes
- `process_gold_standard.py`: Produces the `gold-standard-nodes-uniprot.txt` file from the [EGFR prize file](https://raw.githubusercontent.com/gitter-lab/tps/ca7cafd1e1c17f45ddea07c3fb54d0d70b86ff45/data/resources/eight-egfr-reference-all.txt) from the above paper.
- `map_ensembl.py`: Maps UniProt identifiers to STRING-restricted ENSP identifiers for the STRING-based data.
