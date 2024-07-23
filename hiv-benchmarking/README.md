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


#### 4. filter_empty_pathways.ipynb
- Removes any empty pathways that were created from the SPRAS output in place
- Recommended: make a copy of the directory and then execute code on duplicate directory 


#### 5. plot_num_nodes_summary_table.ipynb
- using the `xx-pathway-summary.txt` file, this notebook produces histograms for the number of nodes present in all the pathways combined
    - i.e. Number of Nodes vs. Count, where count represents the number of total pathways that have a particular number of nodes in their pathway


#### 6. build_kegg-orthology_to_swissprot_map.ipynb
This notebook creates an output of KEGG Orthology mapped to Uniprot and KEGG Orthology mapped to Swissprot, both saved as CSVs. This is done through:
- Using the `biopython` KGMLParser module, the downloaded KGML pathway is parsed to produce a dataframe of proteins (KEGG Orthology IDs)
- Using API calls to the [genome.jp](genome.jp) database, the KEGG orthology IDs are mapped to the KEGG Human protein IDs i.e. the HSA IDs
- Using the `hsa_uniprot.list` file from [LinkDB](https://www.genome.jp/linkdb/) (this was done by downloading the link information between HSA and Uniprot), the HSA IDs are mapped to all the UniProt IDs corresponding to the HSA IDs. Note that this produced a 1-to-many mapping for some proteins
- To remove the 1-to-many mapping and produce a 1-to-1 mapping instead, API calls to [genome.jp](genome.jp) are done again to filter out UniProt IDs that *don't* have SwissProt IDs i.e. protein IDs that haven't been 'manually reviewed'.
- All these mappings are saved locally as csvs.


#### 7. hiv05_comparison_ratios
This notebook contains code to create comparison ratios for SPRAS ensemble pathways vs. the original publication pathway *and* SPRAS ensemble pathways vs. the KEGG pathway 
- This is done by intersecting the nodes found in both pathways, and then comparing the number of nodes found in both to the number of nodes found in the individual original pathways


#### 8. ensemble_node_maxfreq.ipynb
This notebook contains code that loads the ensemble pathway and assigns the highest frequency associated with each node to the respective node.
- The original SPRAS output ensemble pathway file assigns frequencies to edges. 
- The result produced by this notebook assigns the maximum edge frequency associated with each node to the node, creating a max node frequency list that is saved as a CSV.


#### 9. build_prc.ipynb
This notebook produces a precision recall curve using `scikitlearn`. This is done through the following steps:
- Building the PRC uses the max node frequency produced from the ensemble pathway, and the KEGG to Swissprot/Uniprot mapping.
- A vector of 0s and 1s is built that indicates whether the particular protein in the ensemble pathway is found in the KEGG pathway or not. This vector is called `y_kegg` and is attached to the ensemble pathway dataframe.
    - The `max_freq` column in the dataframe is considered as the probabilities/scores and the `y_kegg` column in the dataframe is considered as the true label
- Then, a precision recall curve is built.



