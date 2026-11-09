# spras-met-data
Evaluating potentially suitable metabolic data for Signaling Pathway Reconstruction Analysis Streamliner (SPRAS).
Copied over from irislocal-spras-data, which consisted of:
- Personal tests to evaluate potentially suitable data applications for Signaling Pathway Reconstruction Analysis Streamliner (SPRAS).
    - Specifically regarding application of metabolic signaling pathway reconstruction
- Mainly runs summary statistics through converting PathwayCommons data files from OWL -> SIF -> SPRAS-compatible data formats
- Running tests on SIF formatted files for Recon 2.02, Reactome, KEGG
    - All data is sourced from Pathway Commons at https://www.pathwaycommons.org/pc2/datasources, using v12 files.
        - v12 was used over the latest v14, as the datasets are uploaded as individual files, making them easier to work with.

# Data
## Recon X

Metabolic pathway map of the human body, intended for use as the main background network.

## KEGG

Used as a gold standard network to verify specific metabolic networks to be found within the background network (Recon X). 

# Code
Main files:
- In datasets/kegg/scripts:
    - fetch-from-owl.py
    - kegg-fetchpathway.py
    - listcuratedpathways.py
    - pathway-to-json.py
    - sif-cleaner.py
    - sif-to-spras.py
- In summary-stats:
    - global-sumstats.py
- PathwayCommons has a physical entities file that can be used to progress further on this project (described below), but user will have to download separately due to file size
    - Should be noted that downloading certain files from PathwayCommons such as this one are bulky enough to prevent commits from syncing properly from local -> online repo, use with caution



# Future Directions / Known Issues

- Script needed for recognition of alternate CHEBI ID names for a given compound as the same thing
    - Currently, each small molecule has multiple CHEBI IDs associated with it, preventing meaningful basic analysis on the SPRAS-formatted metabolic datasets
    - Would utilize the physical entities file included in PathwayCommons that lists said compounds + alternate IDs
- Reactome data was initially considered to be implemented as another gold standard metabolic pathway dataset alongside KEGG
    - However, the raw file for it downloaded from PC v12 was far too large to be committed to a repo
    - Also ran into errors when converting it from OWL to SIF format


