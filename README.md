# SPRAS benchmarking

![example workflow](https://github.com/Reed-CompBio/spras-benchmarking/actions/workflows/publish.yml/badge.svg)

Benchmarking datasets for the [SPRAS](https://github.com/Reed-CompBio/spras) project. This repository contains gold standard datasets to evaluate on as well as paper reproductions & improvements incorporating new methodologies.

## Setup

This repository depends on SPRAS. If you want to reproduce the results of benchmarking locally,
you will need to setup SPRAS. SPRAS depends on [Docker](https://www.docker.com/) and [Conda](https://docs.conda.io/projects/conda/en/stable/). If it is hard to install either of these tools,
a [devcontainer](https://containers.dev/) is available for easy setup.

```sh
conda env create -f spras/environment.yml
conda activate spras
pip install ./spras
```

To run the postprocess output scripts, we have a `pyproject.toml` which can be used with your desired python package manager. This separates
the `spras` conda environment from the small scripts we have. (on CI, we use [`uv`](https://docs.astral.sh/uv/).)

To run the benchmarking pipeline, use:

```sh
snakemake --cores 1 --configfile configs/dmmm.yaml --show-failed-logs -s spras/Snakefile
```

> [!NOTE]
> Each one of the dataset categories (at the time of writing, DMMM and PRA) are split into different configuration files.
> Run each one as you would want.
