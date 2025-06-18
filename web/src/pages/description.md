# SPRAS Benchmarking

[SPRAS](https://github.com/Reed-CompBio/spras) is a
utility software designed for performing [signaling pathway](https://en.wikipedia.org/wiki/Cell_signaling#Signal_transduction_pathways) reconstruction
with various algorithms. [SPRAS's documentation](https://spras.readthedocs.io/en/latest/) has more information about its inner workings software.

This benchmarking repository ([see the GitHub](https://github.com/Reed-CompBio/spras-benchmarking/)) is meant to display the performance
of all of the algorithms currently supported by SPRAS on signaling pathways and diseases (to test out <abbr title="Disease Module Mining Methods">DMMMs</abbr>),
comparing their reconstructions with manually curated golden datasets, or synthetically provided datasets from various databases.

All information provided is orchestrated through our GitHub Actions pipeline, and heavy processing is soon to be moved to [HTCondor](https://htcondor.org/).

## Format

Each run's slug has the <code class="color-1">type</code>, the <code class="color-2">dataset</code>, the
<code class="color-3">algorithm</code>, and the <code class="color-4">paramaters</code> [hash](https://en.wikipedia.org/wiki/Hash_function).
