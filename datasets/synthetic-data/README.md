# Synthetic Data

## Download STRING Human Interactome
1. Download the STRING *Homo sapiens* `9606.protein.links.full.v12.0.txt.gz` database file from [STRING](https://string-db.org/cgi/download?sessionId=bL9sRTdIaUEt&species_text=Homo+sapiens&settings_expanded=0&min_download_score=0&filter_redundant_pairs=0&delimiter_type=txt).
2. Move the downloaded file into the `raw/human-interactome/` folder.
3. From the `raw/synthetic-data/` directory, extract the file using:

   ```
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

[Sources](https://www.pnas.org/doi/full/10.1073/pnas.1808790115) are silico human surfaceomes receptors.

[Targets](https://academic.oup.com/nar/article/51/D1/D39/6765312) are human transcription factors. 

## Steps to Generate SPRAS-Compatible Pathways

### 1. Process PANTHER Pathways

1. Open `process_panther_pathway.py` and add the name of any new pathways to the `pathways` vector on **line 6**.
2. Run the command:
   ```
   uv run src/process_panther_pathway.py
   ```
3. This will create five new files in each subfolder of the `pathway-data/` directory:
- `EDGES.txt`
- `NODES.txt`
- `PRIZES-100.txt`
- `SOURCES.txt`
- `TARGETS.txt`

### 2. Convert Pathways to SPRAS-Compatible Format
1.	In `SPRAS_compatible_files.py`, add the name of any new pathways to the `pathway_dirs` list on **line 8**.
2.	From the synthetic-data/ directory, run the command:
```
python src/SPRAS_compatible_files.py
```
3. This will create a new folder named `spras-compatible-pathway-data`, containing subfolders for each PANTHER pathway in SPRAS-compatible format.  
Each subfolder will include the following three files:
- `<pathway_name>_gs_edges.txt`
- `<pathway_name>_gs_nodes.txt`
- `<pathway_name>_node_prizes.txt`

4. From the synthetic-data/ directory, run the command:
```
python src/ratios.py
```
5. This will create a new file `data_ratio.txt` in `spras-compatible-pathway-data` to explain the edge to target/sources ratios.

## Steps to get the interactomes
### 1. Steps to get threshold interactomes
1. From the synthetic-data/ directory, run the command:
```
python src/threshold_interactomes.py
```
2.	This will create a new folder named `interactomes`, containing a subfolder called `uniprot-threshold-interactomes`.
The subfolder will include the following 12 files:
- 10 thresholded interactomes: `uniprot_human_interactome_<threshold>.txt` (thresholds range from 1 to 900)
- `proteins_missing_aliases.csv`: STRING IDs that are missing UniProt accession identifiers
- `removed_edges.txt`: All edges removed from the uniprot_human_interactome_<threshold>.txt files

### 2. Steps to get combined interactomes (Panther pathways and threshold interactomes)
1. In `combine.py`, adjust the `pathway_dirs` list on **line 11** to be the pathways to be included in the combined networks
2. From the synthetic-data/ directory, run the command:
```
python src/combine.py
```
3. This will create a new a subfolder called `uniprot-combined-threshold-interactomes` in `interactomes`.
This subfolder will include 12 files:
- 10 combined threshold interactomes combined with the chosen pathways: `uniprot_combined_interactome_<threshold>.txt` (thresholds range from 1 to 900)
- `overlap_combined_info.csv`
- `overlap_info.csv`

# Pilot Data
For the pilot data, use the list `["Wnt_signaling", "JAK_STAT_signaling", "Interferon_gamma_signaling", "FGF_signaling", "Ras"]` in both:
- the list in `combine.py`
- the list in `overlap_analytics.py`

Make sure these pathways in the list are also added `["Wnt_signaling", "JAK_STAT_signaling", "Interferon_gamma_signaling", "FGF_signaling", "Ras"]`to:
- the `pathways` vector in `ProcessPantherPathway.R`
- the list in `SPRAS_compatible_files.py`

**Once you’ve updated the pathway lists in all relevant scripts, run all the steps above to generate the Pilot dataset.**