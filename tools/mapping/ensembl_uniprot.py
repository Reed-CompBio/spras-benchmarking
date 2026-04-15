import pandas
import os

"""
Utilities for mapping Ensembl and UniProt.

For example,
```py
idmapping_uniprot_ensembl = idmapping_uniprot_mapping(path / "HUMAN_9606_idmapping_selected.tsv")
```

then you can use the `idmapping_as_ensg_uniprot_mapping` or `idmapping_as_ensp_uniprot_mapping` to restrict the mapping to specifically
ENSG or ENSP.
"""


def handle_ensembl_list(idmapping_df: pandas.DataFrame, column_name: str) -> pandas.DataFrame:
    idmapping_df = idmapping_df[idmapping_df[column_name].notnull()]
    # Handle our ;-delimited list
    idmapping_df[column_name] = idmapping_df[column_name].str.split("; ")
    idmapping_df = idmapping_df.explode(column_name)
    # Drop isoforms
    idmapping_df[column_name] = idmapping_df[column_name].str.split(".").str[0]
    idmapping_df = idmapping_df.reset_index(drop=True)
    return idmapping_df


def idmapping_uniprot_mapping(path: str | os.PathLike) -> pandas.DataFrame:
    """
    Gets the UniProt mapping file (`*_idmapping_selected`) as a dataframe with columns
    UniProtKB-AC: High-quality UniProt IDs
    Ensembl: ENSG
    Ensembl_PRO: ENSG (Ensembl Protein IDs)
    """
    # The very powerful UniProt-provided mapping file: its Ensembl mappings are a semicolon-delimited list of Emsembl IDs containing
    # attached isoforms (and not all UniProtKB-AC identifiers have those!) so we'll need to do some extra post-processing.
    # This is `*_idmapping_selected`.
    idmapping_selected_df = pandas.read_csv(
        path,
        header=None,
        # See directory.py for the README associated with this mapping file.
        usecols=[0, 1, 18, 20],
        names=["UniProtKB-AC", "UniProtKB-ID", "Ensembl", "Ensembl_PRO"],
        sep="\t",
    )
    idmapping_selected_df = handle_ensembl_list(idmapping_selected_df, "Ensembl")
    idmapping_selected_df = handle_ensembl_list(idmapping_selected_df, "Ensembl_PRO")
    return idmapping_selected_df


def idmapping_as_ensg_uniprot_mapping(uniprot_mapping: pandas.DataFrame):
    return uniprot_mapping.drop(columns=["Ensembl_PRO"])


def idmapping_as_ensp_uniprot_mapping(uniprot_mapping: pandas.DataFrame):
    return uniprot_mapping.drop(columns=["Ensembl"]).rename(columns={"Ensembl_PRO": "Ensembl"})
