import pandas as pd
import re
import os
from pathlib import Path

"""Extracts gene symbols from dataset for UniProt mapping"""

# configuration - change as needed
input_file = "OmicsSomaticMutationsMatrixDamaging.csv"

dir_path = Path(os.path.dirname(os.path.realpath(__file__)))


def extract_gene_symbols(input_df: pd.DataFrame, input_filename):
    """
    Extracts gene symbols from the input file
    """
    gene_columns = input_df.columns.tolist()[1:]
    gene_symbols = [re.match(r"^(.*?) \(", col).group(1) if " (" in col else col for col in gene_columns]

    gene_symbols_df = pd.DataFrame(gene_symbols, columns=["GeneSymbol"])

    print(f"Extracted {len(gene_symbols_df)} gene symbols from the input file.")
    print(f"First 5 gene symbols: {gene_symbols_df['GeneSymbol'].head().tolist()}")

    return gene_symbols_df


def save_gene_symbols(gene_symbols_df):
    """
    Saves the extracted gene symbols to a CSV file for UniProt web service
    """
    output_dir = dir_path / ".." / "processed"
    os.makedirs(output_dir, exist_ok=True)
    output_path = output_dir / f"DamagingMutationsGeneSymbols.csv"
    gene_symbols_df.to_csv(output_path, index=False)
    return output_path


def main():
    # load dataset
    input_df = pd.read_csv(dir_path / ".." / "raw" / input_file, index_col=0)

    gene_symbols_df = extract_gene_symbols(input_df, input_file)

    output_path = save_gene_symbols(gene_symbols_df)
    print(f"Gene symbols extraction saved to: {output_path}")


if __name__ == "__main__":
    main()
