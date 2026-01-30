# Contributing a new benchmarking dataset

This guide walks new contributors through the process of adding a new dataset for SPRAS benchmarking and running SPRAS on that dataset. It is considered a companion to the [SPRAS Contributing Guide](https://spras.readthedocs.io/en/latest/contributing/index.html), which walks users through adding and algorithm to the SPRAS software. It is useful, but not necessary, to complete that contributing guide before beginning this one.

## Prerequisites

Before following this guide, a contributor will need

- Familiarity with Python ([Carpentries introduction](https://swcarpentry.github.io/python-novice-inflammation/))
- Familiarity with Git and GitHub ([Carpentries introduction](https://swcarpentry.github.io/git-novice/))
- Snakemake ([Carpentries introduction](https://carpentries-incubator.github.io/workflows-snakemake/)
  or [beginner's guide](http://ivory.idyll.org/blog/2023-snakemake-slithering-section-1.html))
- The ability to post files to Google Drive

## Step 0: Fork the repository and create a branch

From the [spras-benchmarking repository](https://github.com/Reed-CompBio/spras-benchmarking),
click the "Fork" button in the upper right corner to create a copy of
the repository in your own GitHub account. Do not change the "Repository
name". Then click the green "Create fork" button.

The simplest way to set up SPRAS benchmarking for local development is to clone your
fork of the repository to your local machine. You can do that with a
graphical development environment or from the command line. After
cloning the repository, create a new git branch called
``example-dataset`` for local neighborhood development. In the
following commands, replace the example username ``agitter`` with your
GitHub username.

```sh
git clone https://github.com/agitter/spras-benchmarking.git
git checkout -b example-dataset
```

Then you can make commits and push them to your fork of the repository
on the ``example-dataset`` branch

```sh
git push origin example-dataset
```

For this example dataset only, you will not merge the changes
back to the original SPRAS benchmarking repository. Instead, you can open a pull
request to your fork so that the SPRAS benchmarking maintainers can still provide
feedback. For example, use the "New pull request" button from
https://github.com/agitter/spras-benchmarking/pulls and set ``agitter/spras-benchmarking`` as both
the base repository and the head repository with ``example-dataset`` as the compare branch.

The [SPRAS Contributing Guide](https://spras.readthedocs.io/en/latest/contributing/index.html) also provides instructions so you can push changes to both the Reed-CompBio version of spras-benchmarking and your fork. 

### Step 1: Install `uv`

Unlike in the main SPRAS repository, we use `uv`, an equivalent to `pip`, for running our dataset pipelines. As we will see later in `1.1`,
we still use Conda for running SPRAS itself.

You can follow `uv`'s installation instructions [on their website](https://docs.astral.sh/uv/getting-started/installation/).

### 1.1: Activate the spras environment and install SPRAS as a submdule.

This repository depends on SPRAS. If you want to reproduce the results of running SPRAS on datasets locally,
you will need to setup SPRAS. SPRAS depends on [Docker](https://www.docker.com/) and [Conda](https://docs.conda.io/projects/conda/en/stable/). If it is hard to install either of these tools,
a [devcontainer](https://containers.dev/) is available for easy setup.

```sh
conda env create -f spras/environment.yml
conda activate spras
pip install ./spras
```

## Step 2: Add a dataset

The goal of a dataset is to take raw data and produce data to be fed to SPRAS. In this guide, we will add a dataset that is provided in `datasets/example`.

### 2.1: Generate an example dataset

Generate a fake dataset by running 

```sh
uv run datasets/example/raw_generation.py
```

The following artifacts will be placed in `dataset/example/`:
- `sources.txt`
- `targets.txt`
- `gold-standard.tsv`
- `interactome.tsv`

### 2.2: Place the example dataset on Google Drive

In more realistic scenarios, the data used in other datasets comes from other sources (whether that's supplementary info in a paper, or out of biological databases like UniProt.) These artifacts can be large, and may occasionally be updated, so we store them in Google Drive for caching and download
them when we want to reconstruct a dataset.

Note that the four artifacts above change every time `raw_generation.py` is run. Upload those artifacts to Google Drive in a folder of your choice. Set the Sharing settings so that _Anyone with the link_ can _View_ the file.

Once shared, copying the URL should look something like:

```
https://drive.google.com/file/d/1Agte0Aezext-8jLhGP4GmaF3tS7gHX-h/view?usp=sharing
```

We always drop the entire `/view?...` suffix, and replace `/file/d/` with `/uc?id=`, which turns the URL to a direct download link, which is internally
downloaded with [gdown](https://github.com/wkentaro/gdown). Those post-processing steps should make the URL now look as so:

```
https://drive.google.com/uc?id=1Agte0Aezext-8jLhGP4GmaF3tS7gHX-h
```

### 2.3: Add the example dataset's location to `cache/directory.py`

Now, add a directive to `cache/directory.py` that specifies the location of the example dataset. This should be added as a new (key,value) pair to the `directory` variable. Since the example dataset doesn't have an online URL, this should use `CacheItem.cache_only`, to indicate that no other online database serves this URL.

Your new directive under the `directory` dictionary should look something as so, with one entry for each of the four artifacts:

```python
...,
"ExampleData": {
    "interactome.tsv": CacheItem.cache_only(
        name="Randomly-generated example data interactome",
        cached="https://drive.google.com/uc?id=..."
    ),
    ...
}
```

Step 3: Set up a workflow to run the example dataset
----------

Now, we need to make these files SPRAS-compatible. To do this, we'll set up a `Snakefile`, which will handle downloading the artifacts from the Google Drive links and running any scripts to reformat the artifacts into SPRAS-compatible formats.

In the example dataset, `sources.txt` and `targets.txt` are already in a SPRAS-ready format, but we need to process `gold-standard.tsv` and `interactome.tsv`.

### 3.1: Write a `Snakefile` to fetch datasets 

Navigate to the `dataset/example` directory and create a `Snakefile` with the top-level directives:

```python
# This provides the `produce_fetch_rules` util to allows us to automatically fetch the Google Drive data.
include: "../../cache/Snakefile"

rule all:
    input:
        # The two files we will be passing to SPRAS
        "raw/sources.txt",
        "raw/targets.txt",
        # The two files we will be processing
        "processed/gold-standard.tsv",
        "processed/interactome.tsv"
```

We'll generate four `fetch` rules, or rules that tell Snakemake to download the data we uploaded to Google Drive earlier.

```python
produce_fetch_rules({
    # The value array is a path into the dictionary from `cache/directory.py`.
    "raw/sources.txt": ["Contributing", "sources.txt"],
    # and so on for targets, gold-standard, and interactome:
    # note that excluding these three stops the Snakemake file from working by design!
    ...
})
```

### 3.2: Write code to put example dataset files in a SPRAS-compatible format

Create two scripts that output SPRAS-ready variants of `raw/gold-standard.tsv` and `raw/interactome.tsv` to `processed/`, consulting
the [SPRAS file format documentation](https://spras.readthedocs.io/en/latest/output.html). You can use any dependencies inside the top-level
`pyproject.toml`, though [pandas](https://pandas.pydata.org/) should suffice, and you can test out your scripts with `uv run <script>`,
an installation requirement from Step 1.


> [!TIP]
> While scripts will usually be run through `Snakemake`, they can also be run standalone through `uv run <script>.py`.
> Users not running scripts from `datasets/<dataset>` will encounter path errors, unless you resolve the file's current directory
> to be its current location through `current_dir = Path(__file__).parent.resolve()`.

### 3.3: Write Snakemake rules to produce SPRAS-compatible files

Once you have your scripts, add rules to the `Snakefile` that consume the raw data and produce your processed data. For example:

```py
rule interactome:
    input:
        "raw/interactome.tsv"
    output:
        "processed/interactome.tsv"
    shell:
        "uv run scripts/process_interactome.py"
```

Once you do the same for `gold-standard.tsv`, your dataset recreation pipeline is ready! This will not run SPRAS itself, but it will allow
your processed dataset files to be reproduced. You can test it with `uv run snakemake --cores 1`.

## Step 4: Add the example dataset to the set of benchmark data

To make sure your dataset is run along with all other datasets when benchmarking is run,
you need to run your new `Snakefile` to `run_snakemake.sh` file in the top-level directory, and add it to the appropiate SPRAS configuration in `configs`.

The example dataset inputs indicate that algorithms designed for pathway reconstruction analysis should be run on this example (as opposed to a disease mining analysis, which would not have sources and targets). Therefore, we will add this dataset to be run when pathway reconstruction analysis (PRA) methods are used. The configuration file for these methods is in `configs/pra.yaml`. 

### Adding to `run_snakemake.sh`

Make sure that your `Snakefile` is run inside the top-level `run_snakemake.sh` file.

### Adding to the SPRAS config

Since this is a pathway problem and not a disease mining problem, we'll mutate `configs/pra.yaml`. Add your dataset and gold standard to the configuration. Since this dataset passes in a mix of raw and processed files, it would be best to make the `data_dir` set to `datasets/example`, then refer to individual files when linking node or edge files in the configuration. Under the `datasets` tag, add lines like this:

```yaml
  - label: exampleDataset
    node_files: ["raw/sources.txt", "raw/targets.txt"]
    edge_files: ["processed/interactome.tsv"]
    data_dir: "datasets/example"
```

To test these, use the `conda` environment from the `spras` submodule to run `snakemake` with SPRAS:

```sh
snakemake --cores 1 --configfile configs/pra.yaml --show-failed-logs -s spras/Snakefile
```

## Making contributions

You can now add your own datasets to the `spras-benchmarking` repo, which will be reviewed by the maintainers. **Check that your data provider isn't already a dataset in `datasets`.** There are some datasets that are able to serve more data, and only use a subset of it: these datasets can be extended for your needs. Code contributions will be licensed using the project's MIT license.

If you wish to contribute to the codebase beyond adding datasets, there are `TODOs` that better enhance the reproducibility and accuracy of datasets or analysis of algorithm outputs, as well as
[open resolvable issues](https://github.com/Reed-CompBio/spras-benchmarking/issues).

If you want to add an algorithm to SPRAS, refer to the [SPRAS repository](https://github.com/Reed-CompBio/SPRAS) instead. If you want to test your new algorithm you PRed to SPRAS, you can swap out the `spras` submodule that this repository uses with your fork of SPRAS.
