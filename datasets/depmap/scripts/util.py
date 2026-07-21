import pandas as pd

def build_hgnc_to_ensp(idmapping_df: pd.DataFrame) -> dict[str, str]:
    """
    Build an HGNC symbol -> ENSP id mapping from the STRING aliases file.
    """

    df_hgnc = (
        idmapping_df.loc[idmapping_df["source"] == "Ensembl_HGNC_symbol"]
        .rename(columns={"#string_protein_id": "ENSP", "alias": "HGNC_symbol"})
        .assign(ENSP=lambda x: x["ENSP"].str.removeprefix("9606."))
    )
    return dict(zip(df_hgnc["HGNC_symbol"], df_hgnc["ENSP"]))
