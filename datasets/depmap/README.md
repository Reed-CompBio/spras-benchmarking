# Cancer Dependency Map Dataset

This folder contains the processed data and the scripts for data analysis and preparation on datasets from The Cancer Dependency Map, an initiative led by the Broad Institute to provide large-scale omics data in identifying cancer dependencies/vulnerabilities.

You can read more about DepMap and the projects included here: https://www.broadinstitute.org/cancer/cancer-dependency-map

## Raw Data
You can visit the DepMap all data downloads portal at: https://depmap.org/portal/data_page/?tab=allData
Download the following datasets under the primary files section of DepMap and move them to a directory named `raw` that you create. The dataset descriptions from the website are also included:

Currently used files:

- `OmicsProfiles.csv`: Omics metadata and ID mapping information for files indexed by Profile ID. This dataset is used for mapping cell line names to DepMap model IDs as a basis for data processing. (file URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=OmicsProfiles.csv)
- `CRISPRGeneDependency.csv`: Gene dependency probability estimates for all models in the integrated gene effect. This dataset is used to identify gold standard genes in each cell line, a dependency probability cutoff of 0.5 is currently used to get the genes with considerable impact on the cell line. (file URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=CRISPRGeneDependency.csv)
- `OmicsSomaticMutationsMatrixDamaging.csv`: Genotyped matrix determining for each cell line whether each gene has at least one damaging mutation. A variant is considered a damaging mutation if LikelyLoF == True. (0 == no mutation; If there is one or more damaging mutations in the same gene for the same cell line, the allele frequencies are summed, and if the sum is greater than 0.95, a value of 2 is assigned and if not, a value of 1 is assigned.). This dataset is used to prepare the input prize file. (file URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=OmicsSomaticMutationsMatrixDamaging.csv)

Future extension files:

- `OmicsExpressionProteinCodingGenesTPMLogp1.csv`: Model-level TPMs derived from Salmon v1.10.0 (Patro et al 2017) Rows: Model IDs Columns: Gene names. (file URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=OmicsExpressionProteinCodingGenesTPMLogp1.csv)
- `OmicsCNGeneWGS.csv`: Gene-level copy number data inferred from WGS data only. Additional copy number datasets are available for download as part of the full DepMap Data Release. (file URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=OmicsCNGeneWGS.csv)


## Scripts
Currently contains:
- `local_cell_line_preprocessing.ipynb`: Jupyter notebook for exploratory data analysis and initial pipeline development. Includes CRISPR dependency analysis with multiple thresholds, visualization of gene dependency distributions, UniProt ID mapping workflow (both gene symbols and gene numbers approaches currently), and step-by-step generation of prize input files and gold standard files for individual cell lines.
- `cell_line_processing.py`: General cell line processing pipeline for generating prize input files and gold standard files converted into Python scripts. Should be reproducible for any cell line name, could be further organized and refined.
- `uniprot_mapping.py`: Gene symbol extraction script for UniProt ID mapping preparation. Parses gene symbols from any DepMap dataset column headers (e.g., "GENE_NAME (12345)" format) and saves them as CSV files ready for input to the UniProt web service. Currently used to extract gene symbols from `OmicsSomaticMutationsMatrixDamaging.csv`, but should be compatible with any omics dataset.


Files used for preparing required files:
- `OmicsProfiles.csv` used for mapping cell line names to DepMap model IDs.
- `OmicsSomaticMutationsMatrixDamaging.csv` used for preparing prize input file.
- `CRISPRGeneDependency.csv` used for preparing gold standard output.  

## Processed Data
Files used for UniProt ID mapping:
- `DamagingMutationsGeneSymbols_20250718.csv`: Gene symbols parsed from gene columns in `OmicsSomaticMutationsMatrixDamaging.csv` on the date described
- `DamagingMutations_idMapping_20250718.tsv`: Gene symbols from `DamagingMutationsGeneSymbols_20250718.csv` mapped to UniProt IDs using UniProt Web Service on the date described
- Folder of processed data for an attempt to do UniProt mapping with the gene index numbers instead, got stuck due to duplicate matches for the same gene number. A future step could be referring to the original mutations file (OmicsSomaticMutations.csv on DepMap, URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=OmicsSomaticMutations.csv) for gene numbers with duplicate matches and do exact matches by seeing where the mutation is located and get more accurate mappings. Contains preliminary processed data (all as of 07/24/2025):
    - `gene_index_mapping_attempt\gene_numbers.txt`: Gene index numbers parsed from gene columns in `OmicsSomaticMutationsMatrixDamaging.csv`
    - `raw_uniprot_idmapping_2025_07_24.tsv`: Initial mapping results, contains both reviewed and unreviewed results, wasn't able to filter directly on UniProt Web Service due to volume
    - `reviewed_id_mapping_2025_07_24.tsv`: Filtered mapping results to only reviewed matches
    - `duplicated_mapping_entries.tsv`: Gene index numbers with duplicate matches

Started processing with the FADU cell line:
- Input prize file prepared from the damaging mutations dataset
- Gold standard file prepared from the CRISPR gene dependency dataset

## config
Example Config file used to get preliminary results on OmicsIntegrator1 and 2 following the EGFR dataset example. Will test out more parameters and update.
The input edge file for the background network can be obtained from the SPRAS repo [`input/phosphosite-irefindex13.0-uniprot.txt`](https://github.com/Reed-CompBio/spras/blob/b5d7a2499afa8eab14c60ce0f99fa7e8a23a2c64/input/phosphosite-irefindex13.0-uniprot.txt)

## Release Citation
For DepMap Release data, including CRISPR Screens, PRISM Drug Screens, Copy Number, Mutation, Expression, and Fusions:
DepMap, Broad (2025). DepMap Public 25Q2. Dataset. depmap.org
