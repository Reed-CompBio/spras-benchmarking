# Contributing

## Helping Out

There are `TODOs` that better enhance the reproducability of datasets or analysis of algorithm outputs, as well as
[open resolvable issues](https://github.com/Reed-CompBio/spras-benchmarking/).

## Adding a dataset

To add a dataset (see `datasets/yeast-osmotic-stress` as an example of a dataset):
1. Check that your dataset provider isn't already added (some of these datasets act as providers for multiple datasets)
1. Create a new folder under `datasets/<your-dataset>`
1. Add a `raw` folder containing your data
1. Add an attached python script that converts your `raw` data to `processed` data
1. If your dataset is a paper reproduction, add a `reproduction/raw` and `reproduction/processed` folder
1. Add your datasets to the appropiate `configs`

## Adding an algorithm

If you want to add an algorithm, refer to the [SPRAS repository](https://github.com/Reed-CompBio/SPRAS) instead.
If you want to test your new algorithm you PRed to SPRAS, you can swap out the `spras` submodule that this repository uses
with your fork of SPRAS.
