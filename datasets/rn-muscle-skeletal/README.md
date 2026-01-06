# ResponseNet Muscle Skeletal dataset.

<!-- TODO: Can we find the source? -->
**This is a paper reproduction**. Other algorithms may also run on this dataset, but no clear source of this has been found.

ResponseNet has not been very clear about the data that they worked on. In this folder, all the files needed to replicate their work is found.

#### File Breakdown
The two files downloaded directly from ResponseNet are `ResponseNetNetwork.json` and `Muscle_Skeletal-Dec2018.tsv`. The JSON file is an output from ResponseNet's sample output, and is what we used to compare to SPRAS.

`sources.txt` and `targets.txt` were manually curated from `ResponseNetNetwork.json`.

The `Muscle_Skeletal-Dec2018.tsv` is the interactome that ResponseNet uses, they do provide a direct download on their site.

#### Other information
In order to download the files for yourself, you can do so at: https://netbio.bgu.ac.il/respnet/, specifically https://netbio.bgu.ac.il/labwebsite/the-responsenet-v-3-web-server-download-page/.

You can directly download the interactome by selecting which one you are interested in using. In order to download their sample, you need to look for the link for the `sample output` and wait for ResponseNet to run. At the time of writing, ResponseNet will not allow you to directly download the source and target files, you must go to the cytoscape section of the software, and download the cytoscape `.json` file.
