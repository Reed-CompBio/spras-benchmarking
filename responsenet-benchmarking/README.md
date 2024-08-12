# ResponseNet Benchmarking Overview

This directory contains all the benchmarking data for the new SPRAS implementation of ResponseNet. 

## Structure of Directory

The sub-directory `/files` are `.json` outputs from the most recent version of the Yeger-Lotem's ResponseNetV.3, and `.txt` outputs from the SPRAS implementation of ResponseNet.

The sub-directory `/output` contains various output files from the `prep_cyto.py` script. For each run of the script, two outputs are created. A nodefile and an edgefile.

Finally, the sub-directory `/images` contains graph exports from cytoscape.

### Original Paper:
Yeger-Lotem E, Riva L, Su LJ, Gitler AD, Cashikar AG, King OD, Auluck PK, Geddie ML, Valastyan JS, Karger DR, Lindquist S, Fraenkel E. Bridging high-throughput genetic and transcriptional data reveals cellular responses to alpha-synuclein toxicity. Nat Genet. 2009 Mar;41(3):316-23. doi: 10.1038/ng.337. Epub 2009 Feb 22. PMID: 19234470; PMCID: PMC2733244.

## Conda Environment

The easiest way to install Python and required packages is with Anaconda.

After installing Anaconda, you can run the following commands in this directory:
```
conda env create -f environment.yml
conda activate responsenet_benchmarking
```
This will create a conda environment that is modified from the SPRAS environment, and then activate the environment. If you want to do it manually, install the latest versions of the following packages: `networkx`, `pybiomart`, and `pandas`.

## Scripts
There is only one script that needs to be run. This is assuming that files have been generated from ResponseNetV.3 and SPRAS. This script is: `prep_cyto.py`. You can call the script by running:
```
python prep_cyto.py --json [relative json filepath] --spras [relative spras filepath] --output [relative output filepath]
```

## Cytoscape Setup
You should open the node and edge files into a cytoscape session. I recommend importing a network from an edge-file, and then adding the nodefile through nodetable import. The cytoscape session file in the directory showcases four networks.
