import pandas as pd
import re
import os
from pathlib import Path

# configuration - change cell line list and dependency cutoff as needed
cell_line_names = ["FADU", "BHY", "SCC4", "INVALID_CELLLINE"]
dependency_threshold = 0.5
require_all_datasets = False  # set to true to require all data types

dir_path = Path(os.path.dirname(os.path.realpath(__file__)))


class CellLineProcessingError(Exception):
    """Custom exception for cell line processing errors."""

    pass


def check_cell_line(cell_line_name, omics_profiles, damaging_mutations_df, CRISPR_dependency, omics_expression, omics_cnv):
    match = omics_profiles[omics_profiles["StrippedCellLineName"].str.lower() == cell_line_name.lower()]

    if match.empty:
        raise CellLineProcessingError(f"Cell line '{cell_line_name}' not found in OmicsProfiles.")

    model_id = match.index[0]
    print(f"Found '{cell_line_name}' cell line, model ID: {model_id}")

    # Check required datasets (mutations + CRISPR) - always raise errors if missing
    if model_id not in damaging_mutations_df.index:
        raise CellLineProcessingError(f"Cell line '{cell_line_name}' not found in required damaging mutations matrix")
    else:
        row_mut = damaging_mutations_df.index.get_loc(model_id)
        print(f"{cell_line_name} cell line (Model ID: {model_id}) present in mutations matrix at row {row_mut}")

    if model_id not in CRISPR_dependency.index:
        raise CellLineProcessingError(f"Cell line '{cell_line_name}' not found in required CRISPR dependency data")
    else:
        row_dep = CRISPR_dependency.index.get_loc(model_id)
        print(f"{cell_line_name} cell line (Model ID: {model_id}) present in CRISPR dependencies at row {row_dep}")

    # check other datasets
    if require_all_datasets:
        if model_id not in omics_cnv.index:
            raise CellLineProcessingError(f"Cell line '{cell_line_name}' not found in required omics CNV data")
        else:
            row_cnv = omics_cnv.index.get_loc(model_id)
            print(f"{cell_line_name} cell line (Model ID: {model_id}) present in omics CNV at row {row_cnv}")

        if model_id not in omics_expression.index:
            raise CellLineProcessingError(f"Cell line '{cell_line_name}' not found in required expression data")
        else:
            row_exp = omics_expression.index.get_loc(model_id)
            print(f"{cell_line_name} cell line (Model ID: {model_id}) present in omics expressions at row {row_exp}")

        print(f"All required datasets contain the cell line '{cell_line_name}' (Model ID: {model_id})")
    else:
        # log availability of optional datasets
        cnv_present = model_id in omics_cnv.index
        expression_present = model_id in omics_expression.index

        if cnv_present:
            row_cnv = omics_cnv.index.get_loc(model_id)
            print(f"{cell_line_name} cell line (Model ID: {model_id}) present in omics CNV at row {row_cnv}")
        else:
            print(f"WARNING: {cell_line_name} not found in omics CNV data (optional dataset)")

        if expression_present:
            row_exp = omics_expression.index.get_loc(model_id)
            print(f"{cell_line_name} cell line (Model ID: {model_id}) present in omics expressions at row {row_exp}")
        else:
            print(f"WARNING: {cell_line_name} not found in expression data (optional dataset)")

        print(f"Required datasets contain the cell line '{cell_line_name}' (Model ID: {model_id})")

    return True, model_id


def generate_prize_files(cell_line_name, model_id, damaging_mutations_df, uniprot_map):
    """Generate prize input file for the cell line based on damaging mutations and uniprot mapping."""
    # Extract mutation data for the cell line
    mutation_row = damaging_mutations_df.loc[model_id]
    # mapping gene symbols to Uniprot IDs
    gene_to_uniprot = dict(zip(uniprot_map["From"], uniprot_map["Entry Name"]))
    rows = []
    for col, score in mutation_row.items():
        # extract gene symbols for mapping
        match = re.match(r"^(.*?) \(", col)
        gene_symbol = match.group(1) if match else col
        # only map for gene symbols in uniprot map
        if gene_symbol in gene_to_uniprot:
            uniprot_id = gene_to_uniprot[gene_symbol]
            rows.append([gene_symbol, uniprot_id, score])
    mapped_prizes_df = pd.DataFrame(rows, columns=["GeneSymbol", "UniprotID", "Prize"])

    # prize input file for all genes
    prizes_input_file = mapped_prizes_df[mapped_prizes_df.columns[1:]].rename(columns={"UniprotID": "NODEID", "Prize": "prize"})
    output_path = dir_path / ".." / "processed" / f"{cell_line_name}_cell_line_prizes.txt"
    prizes_input_file.to_csv(output_path, sep="\t", index=False, header=True)
    print(f"Prize file saved for cell line '{cell_line_name}' at: {output_path}")

    # nonzero prizes input file
    nonzero_prizes_input_file = prizes_input_file[prizes_input_file["prize"] > 0]
    nonzero_output_path = dir_path / ".." / "processed" / f"{cell_line_name}_cell_line_prizes_input_nonzero.txt"
    nonzero_prizes_input_file.to_csv(nonzero_output_path, sep="\t", index=False, header=True)
    print(f"Nonzero prize file saved for cell line '{cell_line_name}' at: {nonzero_output_path}")

    return gene_to_uniprot


def process_single_cell_line(
    cell_line_name: str,
    omics_profiles: pd.DataFrame,
    damaging_mutations_df: pd.DataFrame,
    CRISPR_dependency: pd.DataFrame,
    omics_expression: pd.DataFrame,
    omics_cnv: pd.DataFrame,
    uniprot_map: pd.DataFrame,
):
    """Process a single cell line and generate output files."""
    print(f"\n=== Processing cell line: {cell_line_name} ===")

    is_present, model_id = check_cell_line(cell_line_name, omics_profiles, damaging_mutations_df, CRISPR_dependency, omics_expression, omics_cnv)

    if is_present:
        # make prize input files
        gene_to_uniprot = generate_prize_files(cell_line_name, model_id, damaging_mutations_df, uniprot_map)

        # make gold standard file
        generate_gold_standard(cell_line_name, model_id, CRISPR_dependency, gene_to_uniprot, dependency_threshold)
        print(f"Processing for cell line '{cell_line_name}' completed successfully.")
        return True

def generate_gold_standard(cell_line_name, model_id, CRISPR_dependency, gene_to_uniprot, threshold: float):
    """Generate gold standard file for the cell line based on CRISPR dependency and gene to Uniprot mapping."""
    # map Uniprot IDs to gene symbols in the CRISPR dependency data
    cell_line_dependency = CRISPR_dependency.loc[model_id]
    filtered_dependency = cell_line_dependency[cell_line_dependency > threshold]
    mapped_dependency = []
    for gene, dependency in filtered_dependency.items():
        match = re.match(r"^(.*?) \(", gene)
        gene_symbol = match.group(1) if match else gene
        if gene_symbol in gene_to_uniprot:
            uniprot_id = gene_to_uniprot[gene_symbol]
            mapped_dependency.append([gene_symbol, uniprot_id, dependency])

    mapped_dependency_df = pd.DataFrame(mapped_dependency, columns=["GeneSymbol", "UniprotID", "Dependency"])

    # save mapped dependency as gold standard file
    gold_standard = mapped_dependency_df[mapped_dependency_df.columns[1]]
    threshold_str = str(dependency_threshold).replace(".", "_")
    gold_standard_output_path = dir_path / ".." / "processed" / f"{cell_line_name}_gold_standard_thresh_{threshold_str}.txt"
    gold_standard.to_csv(gold_standard_output_path, sep="\t", index=False, header=False)
    print(f"Gold standard file saved for cell line '{cell_line_name}' at: {gold_standard_output_path}")
    print(f"Threshold: {dependency_threshold} Number of genes in gold standard: {len(gold_standard)}")


def main():
    print(f"Processing cell lines: {cell_line_names}")
    print(f"Require all datasets: {require_all_datasets}")

    try:
        # Load raw dataset files
        print("Loading datasets...")
        base_dir = dir_path / ".." / "raw"
        damaging_mutations_df = pd.read_csv(base_dir / "OmicsSomaticMutationsMatrixDamaging.csv", index_col=0)
        omics_profiles = pd.read_csv(base_dir / "OmicsProfiles.csv", index_col=0)
        omics_expression = pd.read_csv(base_dir / "OmicsExpressionProteinCodingGenesTPMLogp1.csv", index_col=0)
        omics_cnv = pd.read_csv(base_dir / "OmicsCNGeneWGS.csv", index_col=0)
        CRISPR_dependency = pd.read_csv(base_dir / "CRISPRGeneDependency.csv", index_col=0)

        # Load uniprot mapping file form gene symbols
        print("Loading UniProt mapping file...")
        uniprot_map = pd.read_csv(dir_path / ".." / "processed" / "DamagingMutations_idMapping_20250718.tsv", sep="\t")

    except FileNotFoundError as e:
        print(f"ERROR: Required data file not found: {e}")
        raise SystemExit(1)
    except Exception as e:
        print(f"ERROR: Failed to load data files: {e}")
        raise SystemExit(1)

    # Process each cell line
    successful_count = 0
    failed_count = 0
    successful_cell_lines = []
    failed_cell_lines = []

    for cell_line_name in cell_line_names:
        success = process_single_cell_line(
            cell_line_name, omics_profiles, damaging_mutations_df, CRISPR_dependency, omics_expression, omics_cnv, uniprot_map
        )
        if success:
            successful_count += 1
            successful_cell_lines.append(cell_line_name)
        else:
            failed_count += 1
            failed_cell_lines.append(cell_line_name)

    # summary stats
    print("\n=== Processing Summary ===")
    print(f"Total cell lines processed: {len(cell_line_names)}")
    print(f"Successfully processed: {successful_count}")
    print(f"Failed to process: {failed_count}")

    if successful_cell_lines:
        print(f"Successfully processed cell lines: {', '.join(successful_cell_lines)}")

    if failed_cell_lines:
        print(f"Failed cell lines: {', '.join(failed_cell_lines)}")

    # Exit with error code if any cell line processing failed
    if failed_count > 0:
        print("Some cell lines failed to process. Check error messages above.")
        raise SystemExit(1)
    else:
        print("All cell lines processed successfully.")


if __name__ == "__main__":
    main()
