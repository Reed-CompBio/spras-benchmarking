# Cancer Dependency Map Dataset

This folder contains the processed data and the scripts for data analysis and preparation on datasets from The [Cancer Dependency Map](https://depmap.org/portal/), an initiative led by the Broad Institute to provide large-scale omics data in identifying cancer dependencies/vulnerabilities.

You can read more about DepMap and the projects included here: https://www.broadinstitute.org/cancer/cancer-dependency-map

## Raw Data
You can visit the DepMap all data downloads portal at: https://depmap.org/portal/data_page/?tab=allData
Download the following datasets under the primary files section of DepMap and move them to a directory named `raw` that you create. The dataset descriptions from the website are also included:

Currently used files:

- `OmicsProfiles.csv`: Omics metadata and ID mapping information for files indexed by Profile ID. This dataset is used for mapping cell line names to DepMap model IDs as a basis for data processing. (file URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=OmicsProfiles.csv)
- `CRISPRGeneDependency.csv`: Gene dependency probability estimates for all models in the integrated gene effect. This dataset is used to identify gold standard genes in each cell line, a dependency probability cutoff of 0.5 is currently used to get the genes with considerable impact on the cell line. (file URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=CRISPRGeneDependency.csv)
- `OmicsSomaticMutationsMatrixDamaging.csv`: Genotyped matrix determining for each cell line whether each gene has at least one damaging mutation. A variant is considered a damaging mutation if LikelyLoF == True. (0 == no mutation; If there is one or more damaging mutations in the same gene for the same cell line, the allele frequencies are summed, and if the sum is greater than 0.95, a value of 2 is assigned and if not, a value of 1 is assigned.). This dataset is used to prepare the input prize file. (file URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=OmicsSomaticMutationsMatrixDamaging.csv)
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
- `OmicsExpressionProteinCodingGenesTPMLogp1.csv`: Model-level TPMs derived from Salmon v1.10.0 (Patro et al 2017) Rows: Model IDs Columns: Gene names.
- `OmicsCNGeneWGS.csv`: Gene-level copy number data inferred from WGS data only. Additional copy number datasets are available for download as part of the full DepMap Data Release.

## Processed Data
Files used for UniProt ID mapping:
- `DamagingMutationsGeneSymbols.csv`: Gene symbols and Gene IDs parsed from gene columns in `OmicsSomaticMutationsMatrixDamaging.csv` on the date described
- `DamagingMutations_idMapping.tsv`: Gene symbols from `DamagingMutationsGeneSymbols_20250718.csv` mapped to UniProt SwissProt IDs, using Gene ID data
to provide more accurate mappings when possible, since gene symbol -> UniProt mappings are not one-to-one mapping. (TODO: some Gene ID -> UniProt
mappings are also not one-to-one: the accuracy could be improved by identifying the gene via the mutations present in the associated matrix.)

Started processing with the FADU cell line:
- Input prize file prepared from the damaging mutations dataset
- Gold standard file prepared from the CRISPR gene dependency dataset

## Config
Example Config file used to get preliminary results on OmicsIntegrator1 and 2 following the EGFR dataset example. Will test out more parameters and update.

## Data types
The variant (i.e. mutation) scores from DepMap are originally described as
> A variant is considered a damaging mutation if LikelyLoF == True. (0 == no mutation; If there is one or more damaging mutations in the same gene for the same cell line, the allele frequencies are summed, and if the sum is greater than 0.95, a value of 2 is assigned and if not, a value of 1 is assigned.)

Our interpretation of this is
- 0: their variant scoring algorithm assesses this variant to be not functional, that is, not damaging, so we ignore it when assigning prizes
- 1: their scoring algorithm assesses the variant to be functional and it occurs less than all the time (AF < 0.95), which can happen in many ways in cancer (there are multiple clones in the tumor, there is a single clone but one copy of the gene has the mutation and one doesn't); there could also be multiple different mutations in the same gene but they all have frequencies that still sum to less than 0.95
- 2: the gene is predicted to be severely functionally impacted because there are multiple different mutations each assessed to be functional, and cumulatively they have allele frequency approaching 1 or more

We retain the difference in score of 1 or 2 because genes with a score of 2 are more important than genes with a score of 1, per the DepMap scoring scheme.
The DepMap [pipeline documentation](https://storage.googleapis.com/shared-portal-files/Tools/25Q3_Mutation_Pipeline_Documentation.pdf) provides more information about how they derive the original scores, such as VEP for variant effect prediction and gnomAD for allele frequencies.


## Release Citation
For DepMap Release data, including CRISPR Screens, PRISM Drug Screens, Copy Number, Mutation, Expression, and Fusions:
DepMap, Broad (2025). DepMap Public 25Q2. Dataset. depmap.org
