# [SPRAS benchmarking](https://reed-compbio.github.io/spras-benchmarking/)

![example workflow](https://github.com/Reed-CompBio/spras-benchmarking/actions/workflows/publish.yml/badge.svg)

Benchmarking datasets for the [SPRAS](https://github.com/Reed-CompBio/spras) project. This repository contains gold standard datasets to evaluate on as well as paper reproductions & improvements incorporating new methodologies.
The results of every benchmarking run are deployed on GitHub pages. [(See the current web output)](https://reed-compbio.github.io/spras-benchmarking/).

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

## Organization

There are five primary folders in this repository:

```
.
├── configs
├── datasets
├── queue
├── spras
└── web
```

- `configs` contains the YAML files used to configure SPRAS
- `datasets` contains the raw data referenced by `configs`.
- `queue` contains utilities for interacting with [Reed-CompBio/spras-benchmarking-queue](https://github.com/Reed-CompBio/spras-benchmarking-queue)
- `spras` is the cloned submodule of [SPRAS](https://github.com/reed-compbio/spras)
- `web` is an [astro](https://astro.build/) app which generates the `spras-benchmarking` [output](https://reed-compbio.github.io/spras-benchmarking/)

The workflow runs as so:

1. For every dataset, run its inner `Snakefile` with [Snakemake](https://snakemake.readthedocs.io/en/stable/). This is orchestrated
through the top-level [`run_snakemake.sh`](./run_snakemake.sh) shell script.
1. Run each config YAML file in `configs/` with SPRAS.
1. Build the website in `web` with the generated `output` from all of the SPRAS runs, and deploy it on [GitHub Pages](https://pages.github.com/).

For more information on how to add a dataset, see [CONTRIBUTING.md](./CONTRIBUTING.md).
