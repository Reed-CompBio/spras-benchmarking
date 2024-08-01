### Raw Data Sources:


#### 1. 5_min_75%_confidence.sif
- Source: https://github.com/gitter-lab/hiv1-aurkb/blob/main/Results/networks/hiv-networks.cys
- Downloaded the 75% confidence-filtered network at 5 min as a .sif file after opening the file in Cytoscape


#### 2. correctedprize_05.csv and correctedprize_060.csv
- todo


#### 3. david5Min.txt and david60Min.txt
- Source: https://github.com/gitter-lab/hiv1-aurkb/tree/main/Results/enrichment


#### 4. hiv05-networks.cys
- todo


#### 5. hsa_uniprot.list
- Downloaded from https://www.genome.jp/linkdb/ 
- Enter the following in the text boxes under the database link diagram:
    - Click database names to download link information between two databases from *hsa* to *uniprot* in *text* format 


#### idmapping_2024_06_26.tsv & idmapping_active_false_2024_06_26.tsv
- Downloaded from https://www.uniprot.org/ by uploading `prize_list_proteins.txt` (result from `generate_protein_mapping_input.ipynb` notebook) and doing the following mapping:
    - Source database: UniProtKB AC/ID
    - Target database: UniProtKB
- Downloaded the resulting job queue result
- `idmapping_active_false_2024_06_26.tsv` were the obsolete entries found that were seperately downloadeed - these proteins are also part of the `idmapping_2024_06_26.tsv` mapping


#### phosphosite-irefindex13.0-uniprot.txt
- Source: https://github.com/Reed-CompBio/spras/blob/master/input/phosphosite-irefindex13.0-uniprot.txt


#### prize_05.csv and prize_060.csv
- `prize_05.csv` source: https://github.com/gitter-lab/hiv1-aurkb/blob/main/Results/base_analysis/prize_05.csv
- `prize_060.csv` source: https://github.com/gitter-lab/hiv1-aurkb/blob/main/Results/base_analysis/prize_060.csv


#### uniprotkb_go_0000398_AND_taxonomy_id_96_2024_07_25.tsv, 
- todo