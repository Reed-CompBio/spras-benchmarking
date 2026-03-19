# datasets

Datasets contains both the raw data (straight from the study/database), as well as Python scripts and an associated Snakemake file
which take all of the raw data and produce SPRAS-compatible data.

## Prior work

Many of the datasets here have been stripped of their extra post-analysis. Here, we provide commit links to the original work.

- [`hiv`](https://github.com/Reed-CompBio/spras-benchmarking/blob/0293ae4dc0be59502fac06b42cfd9796a4b4413e/hiv-benchmarking)
- [`diseases`](https://github.com/Reed-CompBio/spras-benchmarking/tree/3c0155567dbc43278531b91f9173f6d4f4486dd8/datasets/diseases)
- [`depmap`](https://github.com/Reed-CompBio/spras-benchmarking/tree/b332c0ab53868f111cb89cd4e9f485e8c19aa9e3/datasets/depmap)
- [`yeast_osmotic_stress`](https://github.com/Reed-CompBio/spras-benchmarking/tree/8f69dcdf4a52607347fe3a962b753df396e44cda/yeast_osmotic_stress)
- [`responsenet_muscle_skeletal`](https://github.com/Reed-CompBio/spras-benchmarking/tree/33eacc85cc0749c93f497098c5be00051d3b013a)
- [`egfr`](https://github.com/Reed-CompBio/spras-benchmarking/tree/a80683880a3dc8275ebb6efec16758426bb0f01e)

## `explore` folders

To motivate certain decisions made in-code, such as `synthetic_data`'s PANTHER pathway choices, we provide scripts that use live data
to assist in data curation. These folders can also contain exploratory CLIs for motivating e.g. magic constants.
