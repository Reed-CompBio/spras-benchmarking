# EGFR

This dataset represents protein phosphorylation changes in response to epidermal growth factor (EGF) treatment.

The network includes protein-protein interactions from the (unfortunately now defunct) [iRefIndex](http://irefindex.org/)
and kinase-substrate interactions from [PhosphoSitePlus](http://www.phosphosite.org/).

The processed versions of the provided files can be found in the [Temporal Pathway Synthesizer (TPS)](https://github.com/gitter-lab/tps) repository (and its [associated paper](https://doi.org/10.1016/j.celrep.2018.08.085)). We attach scripts to
reproduce some of the TPS data (when reasonable) to better enable reproducibility.

Specifically, we break this dataset processing down into three phases:
- Collecting raw data
- Reproducing the TPS EGFR data
- Applying the final SPRAS adjustments.

## Reproducing TPS EGFR data

### Interactome

We want to generate the
[phosphosite-irefindex13.0-uniprot.txt](https://github.com/gitter-lab/tps/blob/1d716fb5ae402328a4dd4a43ebe5517bfc67bc31/data/networks/phosphosite-irefindex13.0-uniprot.txt)
file from the TPS repository ([origin commit](https://github.com/koksal/tps/blob/52f9f3da752db8b1be6ada5d7e4216c3984fdba5/data/networks/PhosphoSite_iRefIndex13.0_uniprot_overwrite_pcsf.txt)) ([explanation commit](https://github.com/koksal/tps/pull/4/commits/2219d9570fa5fe85bf47a8bbad8853bf649151f7#diff-c5ad97885a9f1a8caa04e64629373a0484c8c150e9abc68bf74a4ec5c8bdb9c7R30-R32)).

???

### Prizes

We want to generate the [egfr-prizes.txt](https://github.com/gitter-lab/tps/blob/1d716fb5ae402328a4dd4a43ebe5517bfc67bc31/data/pcsf/egfr-prizes.txt) file from the TPS repository. This file is dynamically generated using a [`generate_prizes.sh` script](https://github.com/koksal/tps/blob/bb58d6d89e24dbc39e976a02f1e31387dbe17dfb/pcsf/generate_prizes.sh), which we clean up and embed into this repository (see [scripts/generate_prizes.py](scripts/generate_prizes.py)).

The script depends on three files, all inside the TPS [data/timeseries](https://github.com/koksal/tps/tree/bb58d6d89e24dbc39e976a02f1e31387dbe17dfb/data/timeseries)
folder, which also depend on their own raw files.

`firstfile` and `prevfile`, mapped to `p-values-first.tsv` and `p-values-prev.tsv` are processed from the raw data provided by first supplementary data in the TPS paper. However, since the raw data inside the paper is unlikely
to be updated with the same output format, we trust that the second supplementary data (or the processed data) is correct.

---

`mapfile`, or `peptide-mapping.tsv`, is, as quoted by the README attached to the second supplementary data:
> Obtained by mapping the UniProt accession number (e.g. P00533)
> to the UniProt ID (e.g. EGFR_HUMAN, also known as the UniProt entry name).
> _PSEUDONODE is used to denote a peptide that only maps to obsolete UniProt IDs.

The UniProt accession numbers are provided by the `processed.xlsx` file in the first supplementary dataset. While we just ignored the 'rawer' data for `p-values-first.tsv` and `p-values-prev.tsv`, no data from the second supplementary dataset comes with the associated UniProt accession numbers. TODO: use.

## SPRAS adjustments

These files have been lightly modified for SPRAS by:
- Lowering one edge weight that was greater than 1.
- Removing 182 self-edges.
- Removing a `PSEUDONODE` prize.
- Adding a prize of 10.0 to `EGF_HUMAN`, as the only source is `EGF_HUMAN`.
- Considering all proteins with phosphorylation-based prizes as targets.
- Considering all nodes as active.
- Truncate the prizes down to `Y.YYYYYYYYY` with rounding.

See the `Snakefile` for specific SPRAS post-processing scripts.

## Citation

If you use this dataset, please reference the following publications: (if you only use the background interactome, the latter two publications suffice.)

[Synthesizing Signaling Pathways from Temporal Phosphoproteomic Data](https://doi.org/10.1016/j.celrep.2018.08.085).
Ali Sinan Köksal, Kirsten Beck, Dylan R. Cronin, Aaron McKenna, Nathan D. Camp, Saurabh Srivastava, Matthew E. MacGilvray,
Rastislav Bodík, Alejandro Wolf-Yadlin, Ernest Fraenkel, Jasmin Fisher, Anthony Gitter.
*Cell Reports* 24(13):3607-3618 2018.

[iRefIndex: a consolidated protein interaction database with provenance](https://doi.org/10.1186/1471-2105-9-405).
Sabry Razick, George Magklaras, Ian M Donaldson.
*BMC Bioinformatics* 9(405) 2008.

[PhosphoSitePlus, 2014: mutations, PTMs and recalibrations](https://doi.org/10.1093/nar/gku1267).
Peter V Hornbeck, Bin Zhang, Beth Murray, Jon M Kornhauser, Vaughan Latham, Elzbieta Skrzypek.
*Nucleic Acids Research* 43(D1):D512-520 2015.
