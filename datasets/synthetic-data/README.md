# Synthetic Data

## Download STRING Human Interactome
1. Download the STRING *Homo sapiens* `9606.protein.links.full.v12.0.txt.gz` database file from [STRING](https://string-db.org/cgi/download?sessionId=bL9sRTdIaUEt&species_text=Homo+sapiens&settings_expanded=0&min_download_score=0&filter_redundant_pairs=0&delimiter_type=txt).
2. Move the downloaded file into the `raw/human-interactome/` folder.
3. From the `raw/synthetic-data/` directory, extract the file using:

   ```sh
   gunzip human-interactome/9606.protein.links.full.v12.0.txt.gz
   ```

## Download New PANTHER Pathways
1. Visit [Pathway Commons](https://www.pathwaycommons.org/).
2. Search for the desired pathway (e.g., "signaling") and filter the results by the **PANTHER pathway** data source.  
   Example: [Search for "Signaling" filtered by PANTHER pathway](https://apps.pathwaycommons.org/search?datasource=panther&q=Signaling&type=Pathway)
3. Click on the desired pathway and download the **Extended SIF** version of the pathway.
4. In the `raw/pathway-data/` folder, create a new subfolder named after the pathway you downloaded.
5. Move the downloaded Extended SIF file to this new folder (as a `.txt` file). Rename the file to match the subfolder name exactly.

## Sources and Targets

[Sources](http://wlab.ethz.ch/surfaceome/), or `table_S3_surfaceome.xlsx`, (see [original paper](https://doi.org/10.1073/pnas.1808790115))
are silico human surfaceomes receptors.

[Targets]( https://guolab.wchscu.cn/AnimalTFDB4//#/), or `Homo_sapiens_TF.tsv`, (see [original paper](https://doi.org/10.1093/nar/gkac907))
are human transcription factors.

## Steps to Generate SPRAS-Compatible Pathways

### 1. Process PANTHER Pathways

1. Open `Snakefile` and add the name of any new pathways to the `pathways` entry.
2. Run the command:
   ```sh
   uv run scripts/process_panther_pathway.py <pathway>
   ```
3. This will create five new files in the respective `pathway` subfolder of the `pathway-data/` directory:
- `EDGES.txt`
- `NODES.txt`
- `PRIZES-100.txt`
- `SOURCES.txt`
- `TARGETS.txt`

### 2. Convert Pathways to SPRAS-Compatible Format
1.	In `panther_spras_formatting.py`, add the name of any new pathways to the `pathway_dirs` list on **line 8**.
2.	From the synthetic-data/ directory, run the command:
```
python scripts/panther_spras_formatting.py
```
3. This will create a new folder named `spras-compatible-pathway-data`, containing subfolders for each PANTHER pathway in SPRAS-compatible format.  
Each subfolder will include the following three files:
- `<pathway_name>_gs_edges.txt`
- `<pathway_name>_gs_nodes.txt`
- `<pathway_name>_node_prizes.txt`

# Pilot Data
For the pilot data, use the list `["Wnt_signaling", "JAK_STAT_signaling", "Interferon_gamma_signaling", "FGF_signaling", "Ras"]` in both:
- the list in `combine.py`
- the list in `overlap_analytics.py`

Make sure these pathways in the list are also added `["Wnt_signaling", "JAK_STAT_signaling", "Interferon_gamma_signaling", "FGF_signaling", "Ras"]`to:
- the `pathways` vector in `ProcessPantherPathway.R`
- the list in `panther_spras_formatting.py`

**Once you’ve updated the pathway lists in all relevant scripts, run all the steps above to generate the Pilot dataset.**
