# [SPRAS benchmarking](https://reed-compbio.github.io/spras-benchmarking/)

![example workflow](https://github.com/Reed-CompBio/spras-benchmarking/actions/workflows/publish.yml/badge.svg)

Benchmarking datasets for the [SPRAS](https://github.com/Reed-CompBio/spras) project. This repository contains different dataset collections to evaluate all algorithms in SPRAS.
The results of every benchmarking run are deployed on GitHub pages. [(See the current web output)](https://reed-compbio.github.io/spras-benchmarking/).

## Setup

This repository depends on SPRAS. If you want to reproduce the results locally,
you will need to setup SPRAS. SPRAS depends on [Docker](https://www.docker.com/) and [Conda](https://docs.conda.io/projects/conda/en/stable/). If it is hard to install either of these tools,
a [devcontainer](https://containers.dev/) is available for easy setup.

```sh
conda env create -f spras/environment.yml
conda activate spras
pip install ./spras
```

To run the postprocess output scripts, we have a `pyproject.toml` which can be used with your desired python package manager. This separates
the `spras` conda environment from the small scripts we have. (on CI, we use [`uv`](https://docs.astral.sh/uv/).)

To run the benchmarking pipeline, use (this example is specifically for disease module mining):

```sh
snakemake --cores 1 --configfile configs/scores.yaml --show-failed-logs -s spras/Snakefile
```

To run an individual dataset pipeline, run the respective `Snakefile` in the dataset directory using [uv](https://docs.astral.sh/uv/):

```sh
cd datasets/[dataset]
uv run snakemake --cores 1
```

> [!NOTE]
> Each one of the dataset categories (at the time of writing, scores) are split into different configuration files.
> Run each one as you would want.

## Organization

There are six primary folders in this repository:

```
.
├── cache
├── configs
├── datasets
├── spras
└── web
```

`spras` is the cloned submodule of [SPRAS](https://github.com/reed-compbio/spras), `web` is an
[astro](https://astro.build/) app which generates the `spras-benchmarking` [output](https://reed-compbio.github.io/spras-benchmarking/),
`configs` is the YAML file used to talk to SPRAS, and `datasets` contains the raw data. `cache` is utility for `datasets` which provides a convenient
way to fetch online files for further processing.

The workflow runs as so:

1. For every dataset, run its inner `Snakefile` with [Snakemake](https://snakemake.readthedocs.io/en/stable/). This is orchestrated
through the top-level [`run_snakemake.sh`](./run_snakemake.sh) shell script.
1. Run each config YAML file in `configs/` with SPRAS.
1. Build the website in `web` with the generated `output` from all of the SPRAS runs, and deploy it on [GitHub Pages](https://pages.github.com/).
To see how to build the website, go to its [README](./web/README.md).

For more information on how to add a dataset, see [CONTRIBUTING.md](./CONTRIBUTING.md).
