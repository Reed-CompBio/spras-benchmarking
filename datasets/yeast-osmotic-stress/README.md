# Overview

This project is based on a published case study that studies the osmotic stress response of yeast cells by using proteomic data and Omics Integrator 1 to reconstruct pathways representing the cell response. Our [case study paper](https://doi.org/10.1016/j.celrep.2018.08.085) includes a time-series component and applies Omics Integrator first and then the Temporal Pathway Synthesizer second. Here we only examine the Omics Integrator results of that paper.

The set of files here was used to prepare the input Yeast Proteomic Data and conduct analysis on the output run by Omics Integrator. I swapped the files out as I got different output to analyze and compare. This worked because the number of files I worked with was small.

The major assumption here is that a user will copy the SPRAS repo separately and take the input (the prize1_dummies file and ChasmanNetwork-DirUndir.txt file) and config.yaml files here to run them with SPRAS. Then use the output files from SPRAS as the inputs to the notebooks here. I have included my ensemble file and pathway summary files here in order to run my notebooks as I did.

## Environment

All necessary packages are available at the top-level `pyproject.toml`.

## Scripts

The SPRAS_output folder contains my best SPRAS ensemble output, a single parameter combination output pathway with a Beta parameter of 1.75 exactly, and the pathway summary file for the ensemble file. Copy your files in here to analyze your outputs.

1_Dummy_Node_Add.ipynb - Run 1st: Determines the largest prize value within our input prizes file and adds 3 dummy nodes all assigned with the highest prize to our input file. Outputs a new prizes file with the nodes added. Processes raw prizes file into the prize1_dummies file. Use this prize1_dummies file as your input to SPRAS. Note: I determined that the prizes file already contained 2 of the 5 dummy nodes with prizes, because of this I manually appended the other 3 from the dummy.txt file.

2_Node_Summary_Histo.ipynb - Run 2nd: Takes the pathway-summary file and creates a histogram of the node results that were collected with prizes. Helps begin to understand the outputs.

3_Oi1_Output_Eval.ipynb - Run 3rd: Main analysis file. Takes the best resulting ensemble pathway file, the gold standard nodes, the case study edge results, and the case study edge frequencies as input files. Performs various exploratory data analysis and data prep tasks. Main task is performing set overlap between case study edge results and our results. Creates stats for describing the difference. Includes Venn Diagram visualization code too. One key thing with this file is when trying to analyze the single pathway output file (instead of an ensemble file) you will need to change the path to point inside the folder with the single pathway output file.

File_compare.py - Optional: I used this script to compare two network input files I received to confirm they were in fact the same input I needed for my Omics Integrator input. Can be used to compare any two files passed in as paths. This was specific to my case so if you use the input files here you do not need to run this.

## Future Work

One huge factor in why my results may have been different than the original case study has to do with the lack of a dummy node parameter implemented in the SPRAS version of Omics Integrator 1, which allows a user to pass a file with a list of dummy nodes that the algorithm has to start its reconstructions through. This feature has since been added to SPRAS.

In the case study they ran the tuned parameters with a Beta of 1.75 and r of 0.01 (to add edge noise) and generated 1000 forests. In my case Omics integrator doesn't have a way to run multiple outputs with the same parameter combination in order to ensemble the results and look at edge frequencies. My work around was to use `np.linspace` with a range between 1 and 2 and running 250 - 1000 parameter combinations. The idea being to run parameters as close to 1.75 as possible and compare the outputs.

When I tried to run Cytoscape on anything greater than or equal to 250 combinations, it would hang and then crash with a Java heap space error (see [SPRAS issue](https://github.com/Reed-CompBio/spras/issues/171)). More memory would need to be allocated to potentially fix this.
