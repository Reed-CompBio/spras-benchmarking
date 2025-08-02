import pandas as pd
import re
import os

'''extracts gene symbols from dataset for UniProt mapping'''

# configuration - change as needed
input_file = "OmicsSomaticMutationsMatrixDamaging.csv"
output_date = "20250801"

def extract_gene_symbols(input_df, input_filename):
    """
    Extracts gene symbols from the input file
    """
    gene_columns = input_df.columns.tolist()[1:]
    gene_symbols = [re.match(r"^(.*?) \(", col).group(1) if " (" in col else col for col in gene_columns]

    gene_symbols_df = pd.DataFrame(gene_symbols, columns=["GeneSymbol"])

    print(f"Extracted {len(gene_symbols_df)} gene symbols from the input file.")
    print(f"First 5 gene symbols: {gene_symbols_df['GeneSymbol'].head().tolist()}")

    return gene_symbols_df

def save_gene_symbols(gene_symbols_df, output_date):
    """
    Saves the extracted gene symbols to a CSV file for UniProt web service
    """
    output_dir = os.path.join("..", "processed")
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"DamagingMutationsGeneSymbols_{output_date}_test.csv")
    gene_symbols_df.to_csv(output_path, index=False)
    return output_path

def main():
    # load dataset
    base_dir = os.path.join("..", "raw")
    input_df = pd.read_csv(os.path.join(base_dir, input_file), index_col=0)

    gene_symbols_df = extract_gene_symbols(input_df, input_file)

    output_path = save_gene_symbols(gene_symbols_df, output_date)
    print(f"Gene symbols extraction saved to: {output_path}")


if __name__ == "__main__":
    main()