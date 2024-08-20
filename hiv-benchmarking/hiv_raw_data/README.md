### Raw Data Sources:


#### 1. 5_min_75_confidence.sif
- Source: https://github.com/gitter-lab/hiv1-aurkb/blob/main/Results/networks/hiv-networks.cys
- Downloaded the 75% confidence-filtered network at 5 min as a .sif file after opening the file in Cytoscape and removed the '%' char from the filename to prevent any issues with reading file


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


#### uniprotkb_go_0000398_AND_taxonomy_id_96_2024_07_25.tsv, uniprotkb_go_0045944_AND_taxonomy_id_96_2024_07_29.tsv, uniprotkb_go_0051301_AND_taxonomy_id_96_2024_07_29.tsv, & uniprotkb_go_0098609_AND_taxonomy_id_96_2024_07_29.tsv
- Source: Downloaded from https://www.uniprot.org/ 
- Process:
    - Select Advanced in the main text field
    - Select Gene Ontology in the first blue drop down menu and then select Taxonomy in the second
    - Enter a selected Gene Ontology IDs and then selecting the human taxonomy.
    - Download the tsv with the following columns: Entry, Reviewed, Entry Name, Protein names, Gene Names, Gene Ontology (biological process)
- Selected the following GO IDs: 0000398, 0045944, 0051301, 0098609
    - These were present in `david5Min.txt` and `david6Min.txt`