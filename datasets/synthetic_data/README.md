# Synthetic Data

_Synthetic Data_ is a generic dataset label for a class of synthetic pathways provided by [PathwayCommons](https://www.pathwaycommons.org/).
Currently, we only use [PANTHER](https://pantherdb.org/) pathways from PantherCommons, specifically enumerated in `./pathways.jsonc`.

This entire workflow can also be done with `uv run snakemake --cores 1` inside this directory, as like any other dataset.

## Workflow

The workflow follows these steps in order:

## PANTHER Pathway Fetching

PANTHER pathways are fetched from a singular OWL file containing a bundled collection of all pathways. Since the OWL file that
PathwayCommons provides is over 10gb, we have a separate Snakemake workflow, located under `./panther_pathways`, that trims down the OWL file
to only contain pathways from PANTHER.

Inside `scripts/fetch_pathway.py`, we use this intermediately-generated (and cached!) OWL file to individually generate associated OWL and
SIF files for each pathway.

We have a `./util/parse_pc_pathways.py`, which takes a `pathways.txt` provided by PathwayCommons, and allows us to map the
human-readable pathway names into [identifiers.org](https://identifiers.org/) identifiers, which we later trim down
with our provided list of pathway names in `pathways.jsonc` using `list_curated_pathways.py`.

## SIF Pathway Processing

The scripts `process_panther_pathway.py` and `panther_spras_formatting` convert pathways from the fetching step into ones usable by SPRAS, using
external data:
- [Sources](http://wlab.ethz.ch/surfaceome/), or `table_S3_surfaceome.xlsx`, (see [original paper](https://doi.org/10.1073/pnas.1808790115))
are silico human surfaceomes receptors.
- [Targets]( https://guolab.wchscu.cn/AnimalTFDB4//#/), or `Homo_sapiens_TF.tsv`, (see [original paper](https://doi.org/10.1093/nar/gkac907))
are human transcription factors. We map these to UniProt in `map_transcription_factors.py`.

## Interactome Generation

`interactome.py` uses STRING and UniProt data to produce a UniProt-based interactome.

## Thresholding

Using the interactome and processed pathway files, we threshold pathways. TODO write more about this.
