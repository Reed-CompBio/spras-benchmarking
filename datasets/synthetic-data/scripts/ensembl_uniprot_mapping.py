import pandas
import os

"""
Utilities for getting the UniProt <-> ENSP idmapping through:

idmapping_ensg_df = idmapping_as_ensg_uniprot_mapping(path / "HUMAN_9606_idmapping_selected.tsv")
ensg_to_ensp_df = ensg_to_ensp_mapping(path / "9606.protein.aliases.txt")
idmapping_ensp_df = idmapping_as_ensp_uniprot_mapping(idmapping_ensg_df, ensg_to_ensp_df)
"""

def handle_ensembl_list(
        idmapping_df: pandas.DataFrame,
        column_name: str
) -> pandas.DataFrame:
    idmapping_df = idmapping_df[idmapping_df[column_name].notnull()]
    # Handle our ;-delimited list
    idmapping_df[column_name] = idmapping_df[column_name].str.split("; ")
    idmapping_df = idmapping_df.explode(column_name)
    # Drop isoforms
    idmapping_df[column_name] = idmapping_df[column_name].str.split(".").str[0]
    idmapping_df = idmapping_df.reset_index(drop=True)
    return idmapping_df

def idmapping_uniprot_mapping(
        path: str | os.PathLike
    ) -> pandas.DataFrame:
    # The very powerful UniProt-provided mapping file: its Ensembl mappings are a semicolon-delimeted list of Emsembl IDs containing
    # attached isoforms (and not all UniProtKB-AC identifiers have those!) so we'll need to do some extra post-processing.
    idmapping_selected_df = pandas.read_csv(
        path,
        header=None,
        # See directory.py for the README associated with this mapping file.
        usecols=[0, 18, 20],
        names=["UniProtKB-AC", "Ensembl", "Ensembl_PRO"],
        sep="\t",
    )
    idmapping_selected_df = handle_ensembl_list(idmapping_selected_df, "Ensembl")
    idmapping_selected_df = handle_ensembl_list(idmapping_selected_df, "Ensembl_PRO")
    return idmapping_selected_df

def idmapping_as_ensg_uniprot_mapping(uniprot_mapping: pandas.DataFrame):
    return uniprot_mapping.drop(columns=["Ensembl_PRO"])

def idmapping_as_ensp_uniprot_mapping(uniprot_mapping: pandas.DataFrame):
    return uniprot_mapping.drop(columns=["Ensembl"]).rename(columns={"Ensembl_PRO": "Ensembl"})
