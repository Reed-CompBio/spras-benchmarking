# Synthetic-Benchmarking
This section of the benchmarking analysis is interesting in applying sythentic pathway data from curated databases into SPRAS to see how they perform.

## Setup

Make sure you are in the benchmarking-pipeline directory and you have installed [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)


To setup the environment, do the following command:
```
conda env create -f environment.yml
```

This pipeline is built using the Snakemake workflow. You can configure which pathways you want to include in the `config.yml` file. To run the workflow, do the following command:

```
snakemake --cores 1
```


## Important Files/Directories
- `networks\` is a directory that contains all the methods and the pathway files. They all have the same structure where each algorithm and pathway combination must include the `spras.txt` file with the corresponding `panther.txt` file.
- `scripts\` is a directory that contains all the helper python files in order to generate the auc, heatmap, and stats.
- `Snakefile` this file contains all the rules and workflow requirements.