# WikiPathways

[WikiPathways](https://www.wikipathways.org/) is a crowdsourced wiki of biological pathways.

We use [WikiNetworks](https://doi.org/10.1093/bioinformatics/btab699) to convert signaling pathways of interest into directed graphs.
Due to various flaws, we do clone it locally inside this folder.

As per [Growing DAGs](https://doi.org/10.1101/2022.07.27.501737) and the supplementary section of [PathLinker](https://doi.org/10.1038/npjsba.2016.2),
we use the manually curated transcription factor and receptor lists from
https://doi.org/10.1186/1741-7007-7-50, https://doi.org/10.1016/j.cell.2010.01.044, and https://doi.org/10.1038/nrg2538. See more about how in
`special_nodes.py`.

> NOTE: WikiNetworks is still moderately inaccurate. A further direction could be improving its accuracy.

## Processing

After getting the underlying graph with WikiNetworks, we convert interactions between grouped proteins
(represented by WikiNetworks as a strongly connected subgraph) to their pairwise representation.

We standardize all gene names to Ensembl proteins for use in the StringDB background interactome.

(TODO: as a future direction, these networks would be represented much more nicely as signaling hypergraphs!)
