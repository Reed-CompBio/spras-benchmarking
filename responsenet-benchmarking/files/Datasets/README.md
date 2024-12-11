## Information about Input Files

ResponseNet has not been very clear about the data that they worked on. In this folder, all the files needed to replicate their work is found.

#### File Breakdown
The two files downloaded directly from ResponseNet are `ResponseNetNetwork.json` and `Muscle_Skeletal-Dec2018.tsv`. The Json file is an output from ResponseNet's sample output, and is what we used to compare to SPRAS.

the `sources` and `targets` files are derived from the `ResponseNetNetwork.json`, we went through and scraped all the nodes in the json that were annotated as a source or a target node.

The `Muscle_Skeletal-Dec2018.tsv` is the interactome that ResponseNet uses, they do provide a direct download on their site.

#### Other information
In order to download the files for yourself, you can do so at: http://netbio.bgu.ac.il/respnet

You can directly download the interactome by selecting which one you are interested in using. In order to download their sample, you need to look for the link for the `sample output` and wait for ResponseNet to run. At the time of writing, ResponseNet will not allow you to directly download the source and target files, you must go to the cytoscape section of the software, and download the cytoscape `.json` file.