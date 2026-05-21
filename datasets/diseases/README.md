# GWAS-based Disease Gene Prediction

In this dataset collection, we identify a number of disease-related trait-gene associations from a GWAS database (TIGA). This resources is one of many resources that are integrated into the DISEASES database, which predicts disease-gene associations.

Here, we ask: **how well does GWAS data predict disease-gene associations when they are considered prize nodes within a protein interactome?**

The inputs are the GWAS trait-gene associations (TIGA) and the interactome (STRING-DB). The gold standard dataset is the DISEASES database, which uses other sources of evidence (such as co-occurrence in texts) to establish disease-gene associations.

## DISEASES Database

This dataset comes from the [DISEASES Database](https://diseases.jensenlab.org/About). Relevant papers include:
- Grissa et al., [DISEASES 2.0: a weekly updated database of disease–gene associations from text mining and data integration](https://academic.oup.com/database/article/doi/10.1093/database/baac019/6554833). DATABASE 2022.
- Pletscher-Frankild et al., [DISEASES: Text mining and data integration of disease-gene associations](https://www.sciencedirect.com/science/article/pii/S1046202314003831). Methods, 2015.

Lars Juhl Jensen's lab developed and maintains the STRING database (in collaboration with other groups). The DISEASES Database uses a scheme for scoring disease-gene associations in the same manner as the text-mining scores in STRINGv9.1.

Additionally the DISEASES Database is updated weekly, so this is great for getting relevant information but we should take care in specifying dates/times that the database was accessed.

The DISEASES Database has three channels: text mining, knowledge, and Experiments. **We only consider the Text Mining and Knowledge channels when building the gold standard to avoid overlapping data with the inputs.**

The data can be obtained from [their Downloads page](https://diseases.jensenlab.org/Downloads).

## TIGA

The most recent DISEASES paper (Grissa et al. 2022) integrates a GWAS database called Target Illumination by GWAS Analytics (TIGA), also by the Jensen lab:  
- Yang et al., [TIGA: target illumination GWAS analytics](https://academic.oup.com/bioinformatics/article/37/21/3865/6292081). Bioinformatics 2021.

TIGA calculates confidence scores for gene-trait associations across genome-wide association studies. They include both citation-based and SNP-based measures in their confidence scores (called their mean rank scores); we only take their SNP data for the inputs. These SNPs are weighted by the distance inverse exponential to handle linkage disequilibrium described in their paper (called `N_snpw`). The SNPs themselves are collected from an Ensemble pipeline - TIGA does not do any novel mapping.

The TIGA gene-trait association data can be obtained from [their shiny app page](https://unmtid-shinyapps.net/shiny/tiga/), either through
all gene-trait associations, or through their [archived files](https://unmtid-dbs.net/download/TIGA/).

## Disease Ontology

Finally, we use the Disease Ontology to get from gene-trait associations to gene-disease associations by limiting the traits to diseases. The Disease Ontology data can be obtained from [their Downloads page](https://disease-ontology.org/downloads/).

## Putting it all together

We hashed out this pipeline on the whiteboard in July 2026:

![whiteboard-image](figs/DISEASES-board.jpg)

Here is an updated version after April/May 2026 refactoring:

![workflow-image](figs/DISEASES-workflow.png) 

Briefly the steps are:

**A. Gold Standard Dataset Generation** (`scripts/gold_standard.py`):
- Use the _full_ text mining and knowledge channels from DISEASES.
- For every disease-gene association, get the max value from those two channels (we believe the confidence scores aren't averaged, but that would make sense - we should double-check).
- Remove all disease-gene associations that have a confidence score of less than 4 (retain all w/ scores 4 or 5 out of 5). Call these "high confidence disease-gene pairs."
- Then, remove all disease-gene associations for which there are fewer than 10 high confidence disease-gene pairs for a disease.
- :TODO: should we ensure that all genes are in STRING-DB? If so, should we require that there are 10 high-confidence disease-gene pairs _in the interactome_?

```
Reading data...
 11419274 text mining and 98383 knowledge associations.
 1791 text mining and 97540 knowledge associations with confidenceScore >= 4.

Combining Text Mining and Knowledge Channels...
 1884 associations found in both text mining and knowledge channels. Maximum score retained.
 444 associations from text mining only.
 86628 associations from knowledge only.
 88956 total high confidence disease-gene associations.

Filtering diseases...
 The high-confidence disease-gene pairs correspond to 2200 distinct diseases.
 There are 643 diseases with at least 10 high-confidence disease-gene pairs
 There are 84972 high-confidence disease-gene pairs from the 643 diseases
```

(Note: if you use the filtered datasets, you end with 134 diseases with at least 10 high-confidence disease-gene pairs).

**B. GWAS Dataset Creation** (`scripts/trait_gene_assoc.py`):
- Take the TIGA trait-gene associations and the gold standard Disease Ontology (DO) annotations.
- Map the TIGA traits (in EFO/MONDO/OBA ids) to DOID. Call these "DO-gene associations". There will be snp_w scores for every gene. The mapping is done with [EBI's Ontology xRef Service (OxO)](https://www.ebi.ac.uk/spot/oxo/).
- Ensure all genes are present in the STRING interactome (map the `ENSG` ids to `ENSP` ids). We use the STRING-DB alias mapping to ensure that all TIGA trait-gene associations are in the STRING interactome. (there is a benchmark file for the DISEASES database with STRINGv9.1, but we use the most recent STRING version).

```
Reading TIGA file...
. There are 676837 trait-gene snp scores from TIGA.
.  There are 11978 EFO/MONDO/OBA traits.
Reading gold standard file...
.  There are 643 disease ontology IDs from the gold standard.
Querying https://www.ebi.ac.uk/spot/oxo/api/search?size=500 in batches of 20:
.  [0:20/643]
.  [20:40/643]
.  [40:60/643]
.  [60:80/643]
.  [80:100/643]
.  [100:120/643]
.  [120:140/643]
.  [140:160/643]
.  [160:180/643]
.  [180:200/643]
.  [200:220/643]
.  [220:240/643]
.  [240:260/643]
.  [260:280/643]
.  [280:300/643]
.  [300:320/643]
.  [320:340/643]
.  [340:360/643]
.  [360:380/643]
.  [380:400/643]
.  [400:420/643]
.  [420:440/643]
.  [440:460/643]
.  [460:480/643]
.  [480:500/643]
.  [500:520/643]
.  [520:540/643]
.  [540:560/643]
.  [560:580/643]
.  [580:600/643]
.  [600:620/643]
.  [620:640/643]
.  [640:660/643]
.  342/643 DOIDs are mapped.
.  126/342 of mapped DOIDs are in TIGA.
Mapping TIGA entries to DOID:
.  There are 18422 trait-gene snp scores after mapping to DOID and dropping N/As.
.  There are now 121 DROID traits in the TIGA dataset. If this is fewer than above, then multiple TIGA ids map to the same DOID.
Mapping STRING IDs:
.  There are 18368 DOID-gene snp scores.
```

_Note:_ We discussed a version 2 where we also run DO-gene associations for diseases _not_ in the validation set; that's a later project).

**C. SPRAS Inputs** (`prepare_inputs.py`):
- Identify the diseases that have both TIGA data and also appear in the Gold Standard dataset. 
- Each of the diseases will be a separate node prizes dataset. For each disease, convert the snp_w scores into prizes and make a `node-prizes.txt` file. 
- Each of the diseases will have a validation dataset, comprising of the high confidence diseases-gene pairs from the DISEASES text mining and/or knowledge channels. They have a score (a 4 or a 5), but they are all considered "high confidence" and thus a gene set. 

```
Reading gold standard and trait-gene assoc files:
.  Gold Standard: 84972 disease-gene assocations for 643 diseases.
.  Trait-gene associations (TIGA): 18368 trait-gene associations for 121 diseases.
There are 121 diseases that are in both the gold standard and in TIGA.
Done writing prize and gold standard files for 121 diseases.
```

To make a list of prize files to add to the Snakemake file, use the following command:

```
ls prize_files/* GS_files/* | awk '{print "        \""$1"\","}' 
```
