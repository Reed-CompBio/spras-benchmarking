# HIV dataset

This folder contains data pertaining to following paper: **[HIV-1 virological synapse formation enhances infection spread by dysregulating Aurora Kinase B](https://doi.org/10.1371/journal.ppat.1011492)** - Bruce JW et al. (2023)

The study examines human immune cells responding to viral infection as well as the changes that take place inside the already infected cells, which is the focus here.
The data is from protein abundance and phosphorylation experiments, which will be the input to pathway reconstruction.

**Overarching goal:** Recreate published biological case study on HIV data using SPRAS. This will help in identifying nodes i.e. proteins that are relevant to the disease.

## Raw files

See `raw/README.md`.

## File organization

See `Snakefile` for the way that all of the IO files are connected.

1. `fetch.py` - This grabs the score files from https://doi.org/10.1371/journal.ppat.1011492
1. `prepare.py` - This cleans up the prize files in `raw`; specifically to remove duplicates.
1. `name_mapping.py` - Converts from UniProt KB-ACID to UniProt KB to meet in the middle with `kegg_ortholog.py`. We chose UniProt KB for its generality.
1. `spras_formatting.py` - Formats the input files into a SPRAS-ready format.
1. `kegg_orthology.py` - This is used to generate the KEGG ortholog file for gold standards, but this has yet to be finalized.
