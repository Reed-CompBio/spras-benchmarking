# Overview
Not much to see here....
Actually, the set of scripts here was used to prepare the input Yeast Proteomic Data & conduct analysis on the output run by Omics Integrator. All necessary files required by the scripts are in the Scripts folder. I swapped the files out as I got different output to analyze and compare. This worked because the number of files I worked with was small.


# Scripts
All needed files are in the Required_files folder. You may end up swapping files out in order to analyze different outputs

File Compare - I used this script to compare 2 network input files I received to confirm they were in fact the same input I needed for my Omics Integrator input. Can be used to compare any 2 files passed in as paths. 

Dummy_Node_Add - Determines the largest prize value within our input prizes file and adds 3 dummy nodes all assigned with the highest prize to our input file. Outputs a new prizes file with the nodes added. 

Node_Summary_Histo - Takes the pathway-summary file and creates a histogram of the node results that were collected with prizes.

Oi1_Output_Eval - Main analysis file. Takes the best resulting ensemble pathway file, the KEGG gold standard nodes, the case study edge results and the case study edge frequencies as input files. Performs various EDA & Data prep tasks. Main task is performing set overlap between Case study edge results and our results. Creates stats for describing the difference. Includes Venn Diagram Viz code too. 


