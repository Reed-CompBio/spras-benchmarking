# Cancer Dependency Map Dataset 

This folder contains the processed data and the scripts for data analysis  and preparation on datasets from The Cancer Dependency Map, an initiative lead by the Broad Institute to provide large-scale omics data in identifying cancer dependencies/vulnerabilities. 

You can read more DepMap and the projects included here: https://www.broadinstitute.org/cancer/cancer-dependency-map 

## raw data 
You can visit the DepMap all data downloads portal at: https://depmap.org/portal/data_page/?tab=allData 
Download the following datasets under the primary files section and move them to the raw folder, the dataset descriptions from the website is also included : 


- OmicsProfiles.csv: Omics metadata and ID mapping information for files indexed by Profile ID.This dataset is used for mapping cell line names to DepMap model IDs as a basis for data processing. (file URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=OmicsProfiles.csv) 
- CRISPRGeneDependency.csv: Gene dependency probability estimates for all models in the integrated gene effect.
This dataset is used to identify gold standard genes in each cell line, a dependency probability cutoff of 0.5 is currently used to get the genes with considerable impact on the cell line. (file URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=CRISPRGeneDependency.csv) 
- OmicsCNGeneWGS.csv: Gene-level copy number data inferred from WGS data only.Additional copy number datasets are available for download as part of the full DepMap Data Release.(file URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=OmicsCNGeneWGS.csv) 
- OmicsSomaticMutationsMatrixDamaging.csv: Genotyped matrix determining for each cell line whether each gene has at least one damaging mutation. A variant is considered a damaging mutation if LikelyLoF == True. (0 == no mutation; If there is one or more damaging mutations in the same gene for the same cell line, the allele frequencies are summed, and if the sum is greater than 0.95, a value of 2 is assigned and if not, a value of 1 is assigned.). This dataset is used to prepare the input prize file. (file URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=OmicsSomaticMutationsMatrixDamaging.csv) 
- OmicsExpressionProteinCodingGenesTPMLogp1.csv:Model-level TPMs derived from Salmon v1.10.0 (Patro et al 2017) Rows: Model IDs Columns: Gene names. (file URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=OmicsExpressionProteinCodingGenesTPMLogp1.csv) 

## scripts
Currently only the Jupyter notebook file I used to analyze dependency data and do the data processing locally to get the input prize file and gold standards. Should be reproducible for any cell line name, but is not yet organized or refined for GitHub.  

## processed data 
Files used for Uniprot ID mapping: 
- Gene symbols parsed 
- Gene symbols mapped to Uniprot IDs 
- folder of processed data for an attempt to do UniProt mapping with the gene index numbers instead, got stuck due to duplicate matches for the same gene number, a future step could be referring to the original mutations file(OmicsSomaticMutations.csv on DepMap, URL: https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=OmicsSomaticMutations.csv) for gene numbers with duplicate matches and do exact matches by seeing where the mutation is located and get more accurate mappings. 

Started processing with the FADU cell line: 
- input prize file prepared from the damaging mutations dataset 
- gold standard file prepared from the CRISPR gene dependency dataset

## config 
Example Config file used to get preliminary results on OmicsIntegrator1 and 2 following the EGFR dataset example. Will test out more parameters and update. 

