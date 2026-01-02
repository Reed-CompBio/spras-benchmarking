"""
Grabs the list of known transcription factors and receptors from the following three sources: (where post-filtering of the data here is done:)

- https://doi.org/10.1186/1741-7007-7-50: Mapping the human membrane proteome: a majority of the
human membrane proteins can be classified according to function and evolutionary origin. This paper gives us signaling receptors, stored in
raw/10_1186-data.tsv, which is manually derived from 12915_2009_258_MOESM1_ESM.xls in the supplementary information of the paper,
exported as a .tsv from the Data sub-spreadsheet.
From the PathLinker supplementary information:
> In addition, we manually included three members of the CD3-TCR complex (CD3D, CD3E, and CD3G), which serve as receptors for the
> T Cell Receptor pathway that were not present in the published list.

- https://doi.org/10.1016/j.cell.2010.01.044: An Atlas of Combinatorial Transcriptional Regulation in Mouse and Man. The required data is from
supplementary data Table S1: List of Human and Mouse Transcription Factors Included in This Study, where the human transcription factors
are stored in raw/10_1016-mmc1.tsv (we avoid the mouse transcription factor columns).
As per PathLinker, we take all of the listed TF genes (with appropiate mapping).

- https://doi.org/10.1038/nrg2538: A census of human transcription factors: function, expression and evolution.
The data is from Supplementary Table 2, "List of TF-encoding loci classified as 'a', 'b', 'c', 'x', ...,"
labelled as "Supplementary information S3 (TXT 157 kb)" during retrieval, or raw/st2.tsv
We trim as per PathLinker, "Receptor and TR lists":
> The latter classified TRs as 'a', 'b', 'c', 'x', and
> 'other'. We took only TRs classified as 'a', 'b', or 'other' because TRs in these classes have
> experimental evidence of regulatory function in a mammalian organism or were manually
> curated to be TRs.
"""

import pandas as pd
from pathlib import Path
import os
import io

dir_path = Path(os.path.dirname(os.path.realpath(__file__)))

# The following process functions will return ENSG, not ENSP.
# ENSP mapping will happen later.

def process_repectors():
    """Processes data from https://doi.org/10.1186/1741-7007-7-50"""
    # TODO: the original PathLinker supplementary info section claims to have more (2,124) signaling receptors, depsite
    # this paper only reporting 1,352 signaling receptors,

    df = pd.read_csv(dir_path / ".." / "raw" / "10_1186-data.tsv", sep="\t")
    df = df[df["Main Class"] == "Receptors"]
    df = df[["IPI Accession"]]
    # This paper uses the now-defunct [International Protein Index](https://en.wikipedia.org/wiki/International_Protein_Index)
    print(df)
    pass  # TODO


def process_tf1():
    """Processes data from https://doi.org/10.1016/j.cell.2010.01.044, returning a set of transcription factors from the data source."""
    mmc1 = pd.read_csv(dir_path / ".." / "raw" / "10_1016-mmc1.tsv", sep="\t")
    mmc1 = mmc1[["Entrez gene ID (Human)"]]
    mmc1.columns = ["Entrez"]
    # These are in Entrez genes, so we map them to ENSG via entrez-ensg.tsv
    # TODO: generalize (also used elsewhere)
    entres_df = pd.read_csv(dir_path / ".." / "raw" / "entrez-ensg.tsv", sep="\t", header=None, names=["ENSG", "Entrez"], dtype=str)

    # https://stackoverflow.com/a/39605926/7589775: yikes! We need to mutate these columns to be str outside of the constructor.
    mmc1["Entrez"] = mmc1["Entrez"].astype(str)
    entres_df["Entrez"] = entres_df["Entrez"].astype(str)

    mmc1 = mmc1.merge(entres_df, on="Entrez", how="inner")
    return set(mmc1["ENSG"])


def process_tf2() -> set[str]:
    """Processes data from https://doi.org/10.1038/nrg2538, returning a set of transcription factors from the data source."""
    _description, _delim, st2_text = (dir_path / ".." / "raw" / "st2.tsv").read_text().partition("--Table start--\n")
    st2 = pd.read_csv(io.StringIO(st2_text), sep="\t")
    st2 = st2[["Class", "Ensembl ID"]]
    # As per the top-level description:
    st2 = st2[st2["Class"].isin(["a", "b", "other"])]
    return set(st2["Ensembl ID"])


def main():
    print(process_repectors())
    tfs = process_tf1().union(process_tf2())
    pass  # TODO


if __name__ == "__main__":
    main()
