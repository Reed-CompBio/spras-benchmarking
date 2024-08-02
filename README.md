# Overview
The set of files here was used to prepare the input Yeast Proteomic Data & conduct analysis on the output run by Omics Integrator. All necessary files required by the scripts are in the Scripts folder. I swapped the files out as I got different output to analyze and compare. This worked because the number of files I worked with was small.

The major assumption here is that a user will copy the SPRAS repo separately and take the input (the prize1_dummies file & network1 file) & config files here and run them with SPRAS. Then use the output files from SPRAS as the inputs to the notebooks here. I have included my ensemble file & pathway summary files here in order to run my notebooks as I did. 

# Conda Environment

The easiest way to install Python and the required packages is with [Anaconda](https://www.anaconda.com/download/).

After installing Anaconda, you can run the following commands from the `Scripts` Directory
```
conda env create -f environment.yml
conda activate yeast_env
```
to create a conda environment with the required packages and activate that environment.
If you have a different version of Python already, you can install the specified versions of the required packages in your preferred manner instead of using Anaconda.

# Scripts
The Required_files folder contains my best SPRAS ensemble output, a single parameter combination output pathway with a Beta parameter of 1.75 exactly, and the pathway summary file for the ensemble file. Swap your files in here to analyze your outputs. 

`File_compare.py` - Optional: I used this script to compare 2 network input files I received to confirm they were in fact the same input I needed for my Omics Integrator input. Can be used to compare any 2 files passed in as paths. This was specific to my case so if you use the input files here you do not need to run this. 

1_Dummy_Node_Add - Run 1st: Determines the largest prize value within our input prizes file and adds 3 dummy nodes all assigned with the highest prize to our input file. Outputs a new prizes file with the nodes added. USE THIS prize1_dummies FILE AS YOUR PRIZE INPUT FILE TO SPRAS

2_Node_Summary_Histo - Run 2nd: Takes the pathway-summary file and creates a histogram of the node results that were collected with prizes. Helps begin to understand the outputs. 

3_Oi1_Output_Eval - Run last: Main analysis file. Takes the best resulting ensemble pathway file, the gold standard nodes, the case study edge results, and the case study edge frequencies as input files. Performs various exploratory data analysis & data prep tasks. Main task is performing set overlap between case study edge results and our results. Creates stats for describing the difference. Includes Venn Diagram visualization code too. 


# Future Work

One huge factor in why my results may have been different than my case study has to do with the lack of a dummy node parameter implemented in the SPRAS version of Omics Integrator 1, which allows a user to pass a file with a list of dummy nodes that the algorithm has to start its reconstructions through. 

In the case study they ran the tuned parameters with a Beta of 1.75 and r of 0.01 (to add edge noise) and generated 1000 forests. In my case Omics integrator doesn't have a way to run multiple outputs with the same parameter combination in order to ensemble the results and look at edge frequencies. My work around was to use `np.linspace` with a range between 1 and 2 and running 250 - 1000 parameter combinations. The idea being run parameters as close to 1.75 as possible and compare the outputs. 

When I tried to run Cytoscape on anything greater than or equal to 250 combinations, it would hang and then crash with a Java heap space error (see [SPRAS issue](https://github.com/Reed-CompBio/spras/issues/171). More memory would need to be allocated to potentially fix this. 

