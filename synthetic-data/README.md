# Synthetic Data 
## Steps to Process Panther Pathways
1.	In the pathway-data folder (default: pathway-data/), create a new folder named after the signaling pathway.

    If using a different folder than pathway-data, update the directory variable in ProcessPantherPathway.R.

2.	Inside this new folder, add a .txt file with the same name as the folder.
3.	Open ProcessPantherPathway.R and add to the pathways variable the name of the signaling pathway to the vector.
4.	From the synthetic-data folder, run the script using: Rscript src/ProcessPantherPathway.R

## Steps to Process Panther Pathways
1.	From the synthetic-data folder, run the script using: Rscript src/HumanInteractome.R

## Steps to Make Pathways SPRAS Compatible
1. on line 4, add any new signaling pathway folder name to the list in SPRAS_compatible_files.py

    If using a different folder than pathway-data, update the folder in SPRAS_compatible_files.py on line 5

2. From the synthetic-data folder, run the script using: python src/SPRAS_compatible_files.py

