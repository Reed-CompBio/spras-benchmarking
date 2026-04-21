# EGFR

EGFR dataset. This dataset does a lot less processing for raw files, and is mainly focused on creating the new STRING-based interactome.

## Overview

This produces two datasets: one based on the iRefIndex/PhosphoSite directed interactome of
closed-source origin based off of UniProt identifiers, and another dataset based off of the
more updated though undirected STRING interactome with ENSP identifiers.

## Data

_See `Snakefile` for specific file locations._

The following data is from [_Synthesizing Signaling Pathways from Temporal Phosphoproteomic Data_](https://doi.org/10.1016/j.celrep.2018.08.085):
- Phosphoproteomic prize data (`egfr-prizes.txt`)
- Gold standard nodes (`eight-egfr-reference-all.txt`)
- the `iRefIndex` + `PhosphoSitePlus`-based interactome using UniprotKB IDs.

We also use the [StringDB](https://string-db.org/) (preferably latest, currently v12) human interactome and
[UniProt](https://www.uniprot.org) mapping files which are based in ENSP IDs. See `cache/directory.py` and `./Snakefile` for more info on precisely
where these were fetched.

## Scripts

- `process_prizes.py`: produces a `input-nodes-uniprot.txt` from
[egfr-prizes.txt](https://raw.githubusercontent.com/gitter-lab/tps/refs/heads/master/data/pcsf/egfr-prizes.txt),
trimming psuedonodes, marking all nodes as active, and manually injecting the `EGF_HUMAN`
receptor as a source node (and dummy node for OmicsIntegrator1), while making all other nodes targets and dummy nodes.
- `process_interactome.py`: Produces the STRING `interactome.tsv` file from the STRING links file. Note that
`phosphosite-irefindex13.0-uniprot.txt` is irreproducible, as it has a closed-source origin. It is a directed interactome produced with a
combination of the now archived iRefIndex v13 interactome with extra PhosphoSite-provided nodes.
- `process_gold_standard.py`: Produces the `gold-standard-nodes-uniprot.txt` file from the [EGFR prize file](https://raw.githubusercontent.com/gitter-lab/tps/ca7cafd1e1c17f45ddea07c3fb54d0d70b86ff45/data/resources/eight-egfr-reference-all.txt) from the above paper.
- `map_ensembl.py`: Maps UniProt identifiers to STRING-restricted ENSP identifiers for the STRING-based data.
