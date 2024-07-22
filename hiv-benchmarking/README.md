## HIV benchmarking

This folder contains raw data, processed data, and SPRAS results pertaining to the dataset taken from the following research article: **HIV-1 virological synapse formation enhances infection spread by dysregulating Aurora Kinase B** - Bruce JW, Park E, Magnano C, Horswill M, Richards A, Potts G, et al. (2023)

The data is about human immune cells responding to viral infection and contains data from protein abundance and phosphorylation changes experiments, which will be the input to pathway reconstruction.

**Overarching goal:**
- use SPRAS on a published biological case study on HIV data. This will help in identifying nodes i.e. proteins that are relevant to the disease.


#### 1. compare_prizes_network.ipynb

- This notebook performs data analysis that compares the correctedprize files (5 min and 60 min) from the research article to the original prize files:
    - Some proteins in the original prize files have the syntax `majorIdentifier-N ` where N denotes isoforms. 
    - Data analysis involves checking how many proteins in the original prize files are repeats of the same majorIdentifier
    - Additionally, it involves checking if the network file used (`phosphosite-irefindex13.0-uniprot.txt`) contains secondary identifiers with the -N syntax. 
    - If it does not, we will want to strip that syntax from the prize file as part of preprocessing.
