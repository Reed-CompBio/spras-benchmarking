# Synthetic Data

_Synthetic Data_ is a generic dataset label for a class of synthetic pathways provided by [PathwayCommons](https://www.pathwaycommons.org/).
Currently, we only use [PANTHER](https://pantherdb.org/) pathways from PantherCommons, specifically enumerated in `./pathways.jsonc`.

This entire workflow can also be done with `uv run snakemake --cores 1` inside this directory, as like any other dataset.

## Workflow

## PANTHER Pathway Fetching

PANTHER pathways are fetched from a singular OWL file containing a bundled collection of all pathways. Since the OWL file that
PathwayCommons provides is over 10gb, we have a separate Snakemake workflow, located nuder `./panther_pathways`, that trims down the OWL file
to only contain pathways from PANTHER.

Inside `scripts/fetch_pathway.py`, we use this intermediately-generated (and cached!) OWL file to individually generate associated OWL and
SIF files for each pathway.

We have a `./util/parse_pc_pathways.py`, which takes a `pathways.txt` provided by PathwayCommons, and allows us to map the
human-readable pathway names in `pathways.jsonc` into [identifiers.org](https://identifiers.org/) identifiers.

## Sources and Targets

[Sources](http://wlab.ethz.ch/surfaceome/), or `table_S3_surfaceome.xlsx`, (see [original paper](https://doi.org/10.1073/pnas.1808790115))
are silico human surfaceomes receptors.

[Targets]( https://guolab.wchscu.cn/AnimalTFDB4//#/), or `Homo_sapiens_TF.tsv`, (see [original paper](https://doi.org/10.1093/nar/gkac907))
are human transcription factors.

### 1. Process PANTHER Pathways

1. Open `Snakefile` and add the name of any new pathways to the `pathways` entry.
2. Run the command:
   ```sh
   uv run scripts/process_panther_pathway.py <pathway>
   ```
3. This will create five new files in the respective `pathway` subfolder of the `pathway-data/` directory:
- `edges.txt`
- `nodes.txt`
- `prizes-100.txt`
- `sources.txt`
- `targets.txt`

### 2. Convert Pathways to SPRAS-Compatible Format
1.	In `panther_spras_formatting.py`, add the name of any new pathways to the `pathway_dirs` list on **line 8**.
2.	From the synthetic_data/ directory, run the command:
```
python scripts/panther_spras_formatting.py
```
3. This will create a new folder named `spras-compatible-pathway-data`, containing subfolders for each PANTHER pathway in SPRAS-compatible format.  
Each subfolder will include the following three files:
- `<pathway_name>_gs_edges.txt`
- `<pathway_name>_gs_nodes.txt`
- `<pathway_name>_node_prizes.txt`

