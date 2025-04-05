# spras-benchmarking
Benchmarking datasets for the [SPRAS](https://github.com/Reed-CompBio/spras) project

# Synthetic Data 
## Steps to Process Panther Pathways
1.	In the pathway-data directory (default: pathway-data/), create a new folder named after the signaling pathway.

    If using a different directory, update the directory in ProcessPantherPathway.R accordingly.

2.	Inside this new folder, add a .txt file with the same name as the folder.
3.	Open ProcessPantherPathway.R and set the pathway variable to the name of the signaling pathway.
4.	From the synthetic-data directory, run the script using: Rscript src/ProcessPantherPathway.R

## Steps to Process Panther Pathways
1.	From the synthetic-data directory, run the script using: Rscript src/HumanInteractome.R
