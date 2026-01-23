# HIV dataset

This folder contains data pertaining to following paper: **[HIV-1 virological synapse formation enhances infection spread by dysregulating Aurora Kinase B](https://doi.org/10.1371/journal.ppat.1011492)** - Bruce JW et al. (2023)

The study examines human immune cells responding to viral infection as well as the changes that take place inside the already infected cells, which is the focus here.
The data is from protein abundance and phosphorylation experiments, which will be the input to pathway reconstruction.

**Overarching goal:** Recreate published biological case study on HIV data using SPRAS. This will help in identifying nodes i.e. proteins that are relevant to the disease.

## Raw files

Follow the `Snakemake` directive to find the fetched URLs for these.

- `prize_05.tsv`: Prizes files generated after 5 minutes, from the original paper above.
- `prize_060.tsv`: Prizes files generated after 60 minutes, from the original paper above.
- `ko03250.xml`: KEGG Orthology Pathway ID 03250
- `HUMAN_9606_idmapping.tsv`: File provided by UniProt, used for mapping UniProt identifiers for  `name_mapping.py`.
- `phosphosite-irefindex13.0-uniprot.txt`: The background interactome from the now-gone iRefIndex.

## File organization

See `Snakefile` for the way that all of the IO files are connected.

1. `prepare.py` - This cleans up the prize files in `raw`; specifically to remove duplicates, and to prepare the list of UniProt nodes to be mapped by `name_mapping.py`
1. `name_mapping.py` - Converts from UniProt KB-ACID to UniProt KB to meet in the middle with `kegg_ortholog.py`, and to match with the proteins for the iRefIndex interactome. We chose UniProt KB for its generality. We also remove identifiers with an `-N` suffix and remove duplicates, to make sure isoforms aren't considered as distinct during pathway reconstruction.
1. `spras_formatting.py` - Formats the input files into the universal SPRAS format.
1. `kegg_orthology.py` - This is used to generate the KEGG ortholog file for gold standards, but this has yet to be finalized.
