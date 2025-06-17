# spras-benchmarking

Benchmarking datasets for the [SPRAS](https://github.com/Reed-CompBio/spras) project

## Setup

This repository depends on SPRAS. If you want to reproduce the results of benchmarking locally,
you will need to setup SPRAS. SPRAS depends on Docker and Conda - if it is hard to install either,
a devcontainer is available for easy setup.

```sh
conda env create -f spras/environment.yml
conda activate spras
pip install ./spras
```

To run the benchmarking pipeline, use:

```sh
snakemake --cores 1 --configfile configs/config.yaml --show-failed-logs -s spras/Snakefile
```
