# WikiPathways

[WikiPathways](https://www.wikipathways.org/) is a crowdsourced wiki of biological pathways.

We use [WikiNetworks](https://doi.org/10.1093/bioinformatics/btab699) to convert signaling pathways of interest into directed graphs.
Due to various flaws, we do clone it locally inside this folder.

## Processing

After getting the underlying graph with WikiNetworks, we convert interactions between grouped proteins
(represented by WikiNetworks as a strongly connected subgraph) to their pairwise representation.

We standardize all gene names to Ensembl proteins for use in the StringDB background interactome.

(TODO: as a future direction, these networks would be represented much more nicely as signaling hypergraphs!)
