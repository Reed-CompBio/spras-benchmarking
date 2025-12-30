import pandas as pd
import os
from pathlib import Path

"""Extracts gene symbols from dataset for UniProt mapping"""

dir_path = Path(os.path.dirname(os.path.realpath(__file__)))


def extract_gene_symbols(input_df: pd.DataFrame) -> pd.DataFrame:
    """
    Extracts gene symbols from the input file.
    """
    # The gene symbols are stored in the columns of the matrix, as GENE_NAME (GENE_ID)
    gene_columns = input_df.columns.tolist()[1:]
    gene_symbols = [
        # We want to extract GENE_NAME from GENE_NAME (Unknown)
        (col[:col.find("(") - 1], None) if "(Unknown)" in col else
        # or GENE_ID from "GENE_NAME (GENE_ID)"
        (col[:col.find("(") - 1], col[col.find("(") + 1:-1]) if "(" in col else
        (col, None)
        for col in gene_columns
    ]

    gene_symbols_df = pd.DataFrame(gene_symbols, columns=["GeneSymbol", "GeneID"])

    print(f"Extracted {len(gene_symbols_df)} gene symbols from the input file.")
    print(f"First 5 gene symbols: {gene_symbols_df['GeneSymbol'].head().tolist()}")

    return gene_symbols_df


def main():
    # Load the dataset
    # We only read the first row since we only care about the column names of the matrix
    input_df = pd.read_csv(dir_path / ".." / "raw" / "OmicsSomaticMutationsMatrixDamaging.csv", index_col=0, nrows=1)

    gene_symbols_df = extract_gene_symbols(input_df)

    processed_path = dir_path / ".." / "processed"
    processed_path.mkdir(exist_ok=True)

    # We split gene_symbols_df into where GeneID does and does not exist
    gene_symbols_df_id = gene_symbols_df[gene_symbols_df["GeneID"].notnull()]
    gene_symbols_df_nid = gene_symbols_df[~gene_symbols_df["GeneID"].notnull()]

    # Now, we map it through HUMAN_9606_idmapping_selected.tsv from UniProt (see fetch.py):
    # Associated documentation of this format is at https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
    # We have two files: idmapping_selected will be used for GeneID -> UniProtKB-AC mapping,
    # while idmapping will be used for GeneSymbol -> UniProtKB-AC mapping.

    # We'll also take the idmapping data and trim for specifically Swiss-Prot (curated) genes.
    curated_df = pd.read_csv(dir_path / ".." / "raw" / "SwissProt_9606.tsv", sep='\t', usecols=["Entry", "Entry Name", "Gene Names"])
    curated_df.columns = ["UniProtKB-AC", "Entry Name", "Gene Names"]

    idmapping_df = pd.read_csv(dir_path / ".." / "raw" / "HUMAN_9606_idmapping.tsv", header=None, names=["UniProtKB-AC", "ID_type", "Value"], sep='\t')
    idmapping_df = idmapping_df[idmapping_df["ID_type"] == "Gene_Name"].drop(columns=["ID_type"]).rename(columns={"Value": "GeneSymbol"})
    idmapping_df = idmapping_df.merge(curated_df, on="UniProtKB-AC", how="inner")
    gene_symbols_df_nid = gene_symbols_df_nid.merge(idmapping_df, on="GeneSymbol", how="inner").drop(columns=["GeneID"])

    idmapping_selected_df = pd.read_csv(
        dir_path / ".." / "raw" / "HUMAN_9606_idmapping_selected.tsv",
        header=None, usecols=[0, 1, 2], names=["UniProtKB-AC", "UniProtKB-ID", "GeneID"], sep='\t'
    )
    idmapping_selected_df = idmapping_selected_df[~idmapping_selected_df["GeneID"].isna()]
    idmapping_selected_df = idmapping_selected_df.merge(curated_df, on="UniProtKB-AC", how="inner")
    gene_symbols_df_id = gene_symbols_df_id.merge(idmapping_selected_df, on="GeneID", how="inner")
    gene_symbols_df_id = gene_symbols_df_id.drop(columns=["GeneID", "UniProtKB-ID"])

    gene_symbol_df = gene_symbols_df_id.merge(gene_symbols_df_nid, on=["GeneSymbol", "UniProtKB-AC", "Entry Name", "Gene Names"], how="outer")
    gene_symbol_df = gene_symbol_df.drop(columns=["Gene Names"])
    gene_symbol_df = gene_symbol_df.rename(columns={"GeneSymbol": "From"})
    gene_symbol_df.to_csv(dir_path / ".." / "processed" / "DamagingMutations_idMapping.tsv", sep='\t', index=False)


if __name__ == "__main__":
    main()
