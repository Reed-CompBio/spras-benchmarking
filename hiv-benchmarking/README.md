## HIV benchmarking

This folder contains raw data, processed data, and SPRAS results pertaining to the dataset taken from the following research article: **HIV-1 virological synapse formation enhances infection spread by dysregulating Aurora Kinase B** - Bruce JW, Park E, Magnano C, Horswill M, Richards A, Potts G, et al. (2023)

The data is about human immune cells responding to viral infection and contains data from protein abundance and phosphorylation changes experiments, which will be the input to pathway reconstruction.

**Overarching goal:** use SPRAS on a published biological case study on HIV data. This will help in identifying nodes i.e. proteins that are relevant to the disease.


#### 1. compare_prizes_network.ipynb

This notebook performs data analysis that compares the correctedprize files (5 min and 60 min) from the research article to the original prize files:
- Some proteins in the original prize files have the syntax `majorIdentifier-N ` where N denotes isoforms. 
- Data analysis involves checking how many proteins in the original prize files are repeats of the same majorIdentifier
- Additionally, it involves checking if the network file used (`phosphosite-irefindex13.0-uniprot.txt`) contains secondary identifiers with the -N syntax. 
- If it does not, we will want to strip that syntax from the prize file as part of preprocessing.


#### 2. generate_protein_mapping_input.ipynb
This notebook creates the list of proteins to upload to UniProt as a .txt file. 
- The proteins from the the `prize05.csv` and `prize060.csv` files are combined together and this is saved as a csv. The file is saved as `prize_list_proteins.txt`.


#### 3. preprocess_prize_file.ipynb

This notebook preprocesses the original prize files into a SPRAS compatible format. The following preprocessing steps are done:
- Node Identifier Simplification:
    - Original Format: Columns - UniprotID (i.e. protein IDs) and prize. Some proteins had the syntax majorIdentifier-N (N denotes isoforms).
    - Modification: Remove the -N suffix to retain only the majorIdentifier.
    - Prize Selection: Retain the maximum prize number associated with each major protein identifier.
- Protein Mapping:
    - Network file used: `hiv_raw_data/phosphosite-irefindex13.0-uniprot.txt`
    - Issue: Node file and network file use different protein codes.
    - Solution: Map and replace protein identifiers in the node file to match the network file using UniProt database. The mapping file used was `hiv_raw_data/idmapping_2024_06_26.tsv`, which was downloaded from UniProt. To get the most recent UniProt mapping, run `generate_protein_mapping_input.ipynb`, upload that to Uniprot, and then use the resulting `idmapping`.
- Column Header Update:
    - Original Headers: UniProtID and prize
    - New Headers: NODEID and prize
- These modified files are saved as `hiv_processed_data/modified_prize_xx.csv`


#### 




