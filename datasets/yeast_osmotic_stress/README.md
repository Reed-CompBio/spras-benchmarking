# Overview

This dataset is based on [_Synthesizing Signaling Pathways from Temporal Phosphoproteomic Data_](https://doi.org/10.1016/j.celrep.2018.08.085), which studies the osmotic stress response of yeast cells by using proteomic data and Omics Integrator 1 to reconstruct pathways representing the cell response.

**NOTE**: The original paper also included a time-series component to use the [Temporal Pathway Synthesizer](https://doi.org/10.1016/j.celrep.2018.08.085). Here, until SPRAS supports temporal graphs, we only examine the non-temporal parts of this paper - specifically, we aim to reproduce the [OmicsIntegrator1](https://github.com/fraenkel-lab/omicsIntegrator) results.

The set of files here was used to prepare the input Yeast Proteomic Data and conduct analysis on the output run by Omics Integrator. I swapped the files out as I got different output to analyze and compare. This worked because the number of files I worked with was small.

The major assumption here is that a user will copy the SPRAS repo separately and take the input (the prize1_dummies file and ChasmanNetwork-DirUndir.txt file) and config.yaml files here to run them with SPRAS. Then use the output files from SPRAS as the inputs to the notebooks here. I have included my ensemble file and pathway summary files here in order to run my notebooks as I did.

## Scripts

There is only one script, `process_prizes.py`, which:
1. Determines the largest prize value within our input prizes file and adds 3 dummy nodes all assigned with (the highest prize? TODO: This seems to currently a magic value.)
1. Outputs a new prizes file with the nodes added.
1. Processes raw prizes file into the `prize1_dummies.txt` file as SPRAS input. Note: I determined that the prizes file already contained 2 of the 5 dummy nodes with prizes, because of this I manually appended the other 3 from the dummy.txt file.

## Raw Files

There are other raw files inside the Snakefile, but we don't use them here. We focus on these two raw files instead:
- `prizes.txt`: From supplementary data 3, containing the prize data to be fed into SPRAS for reconstruction.
- `ChasmanNetwork-DirUndir.txt`: The background interactome provided by [Pathway connectivity and signaling coordination in the yeast stress‐activated signaling network](https://doi.org/10.15252/msb.20145120).

## Future Work

(_Note: results are from [this `config.yaml`](https://github.com/tristan-f-r/spras-benchmarking/blob/9477d85871024a5e3a4b0b8b9be7e78c0d0ee961/yeast_osmotic_stress/config.yaml)_).
One huge factor in why my results may have been different than the original case study has to do with the lack of a dummy node parameter implemented in the SPRAS version of Omics Integrator 1, which allows a user to pass a file with a list of dummy nodes that the algorithm has to start its reconstructions through. This feature has since been added to SPRAS.

In the case study they ran the tuned parameters with a Beta of 1.75 and r of 0.01 (to add edge noise) and generated 1000 forests. In my case Omics integrator doesn't have a way to run multiple outputs with the same parameter combination in order to ensemble the results and look at edge frequencies. My work around was to use `np.linspace` with a range between 1 and 2 and running 250 - 1000 parameter combinations. The idea being to run parameters as close to 1.75 as possible and compare the outputs.

When I tried to run Cytoscape on anything greater than or equal to 250 combinations, it would hang and then crash with a Java heap space error (see [SPRAS issue](https://github.com/Reed-CompBio/spras/issues/171)). More memory would need to be allocated to potentially fix this.
