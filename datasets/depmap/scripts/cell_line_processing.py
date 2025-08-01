import pandas as pd
import re
import os

# configuration - change cell line and dependency cutoff as needed
cell_line_name = "FADU"
dependency_threshold = 0.5


def check_cell_line(cell_line_name, omics_profiles, damaging_mutations_df, CRISPR_dependency, omics_expression, omics_cnv):
    print(f"Checking cell line: {cell_line_name} presence in all datasets...")
    match = omics_profiles[omics_profiles["StrippedCellLineName"].str.lower()== cell_line_name.lower()]


    if match.empty:
        print(f"Cell line '{cell_line_name}' not found in OmicsProfiles.")
        return False, None

    model_id = match.index[0]
    profile_row = match.index[0]
    print(f"Found '{cell_line_name}' cell line, model ID: {model_id} (row {profile_row})")

    # check other datasets
    dataset_check = []
    # check with damaging mutations matrix
    if model_id in damaging_mutations_df.index:
        row_mut = damaging_mutations_df.index.get_loc(model_id)
        print(f"{cell_line_name} cell line (Model ID: {model_id}) present in mutations matrix at row {row_mut} ")
        dataset_check.append(True)
    else:
        print(f"{cell_line_name} not found in mutation matrix")
        dataset_check.append(False)
    # check with CRISPR dependency
    if model_id in CRISPR_dependency.index:
        row_dep = CRISPR_dependency.index.get_loc(model_id)
        print(f"{cell_line_name} cell line (Model ID: {model_id}) present in CRISPR dependencies at row {row_dep} ")
        dataset_check.append(True)
    else:
        print(f"{cell_line_name} not found in CRISPR dependency data")
        dataset_check.append(False)
    # check with CNV
    if model_id in omics_cnv.index:
        row_cnv = omics_cnv.index.get_loc(model_id)
        print(f"{cell_line_name} cell line (Model ID: {model_id}) present in omics CNV at row {row_cnv} ")
        dataset_check.append(True)
    else:
        print(f"{cell_line_name} not found in omics CNV data")
        dataset_check.append(False)
    #check with expression
    if model_id in omics_expression.index:
        row_exp = omics_expression.index.get_loc(model_id)
        print(f"{cell_line_name} cell line (Model ID: {model_id}) present in omics expressions at row {row_exp} ")
        dataset_check.append(True)
    else:
        print(f"{cell_line_name} not found in expression data ")
        dataset_check.append(False)

    all_present = all(dataset_check)
    if all_present:
        print(f"All datasets contain the cell line '{cell_line_name}' (Model ID: {model_id})")
    else:
        print(f"{cell_line_name} cell line is missing in some datasets. Cannot proceed. Check the output above for details.")
    return all_present, model_id

def generate_prize_files(cell_line_name, model_id, damaging_mutations_df, uniprot_map):
    '''Generate prize input file for the cell line based on damaging mutations and uniprot mapping.'''
    print(f"Generating prize input file for cell line '{cell_line_name}'...")
    # Extract mutation data for the cell line
    mutation_row = damaging_mutations_df.loc[model_id]
    # mapping gene symbols to Uniprot IDs
    gene_to_uniprot = dict(zip(uniprot_map["From"], uniprot_map["Entry Name"]))
    rows = []
    for col, score in mutation_row.items():
        #extract gene symbols for mapping
        match = re.match(r"^(.*?) \(", col)
        gene_symbol = match.group(1) if match else col
        # only map for gene symbols in uniprot map
        if gene_symbol in gene_to_uniprot:
            uniprot_id = gene_to_uniprot[gene_symbol]
            rows.append([gene_symbol, uniprot_id, score])
    mapped_prizes_df = pd.DataFrame(rows, columns=["GeneSymbol", "UniprotID", "Prize"])

    # prize input file for all genes
    prizes_input_file = mapped_prizes_df[mapped_prizes_df.columns[1:]].rename(columns={"UniprotID": "NODEID", "Prize":"prize"})
    output_path = os.path.join("..", "processed", f"{cell_line_name}_cell_line_prizes_all.txt")
    prizes_input_file.to_csv(output_path, sep='\t', index=False, header=True)
    print(f"Prize file saved for cell line '{cell_line_name}' at: {output_path}")

    # nonzero prizes input file
    nonzero_prizes_input_file = prizes_input_file[prizes_input_file["prize"] > 0]
    nonzero_prizes_input_file
    nonzero_output_path = os.path.join("..", "processed", f"{cell_line_name}_cell_line_prizes_input_nonzero.txt")
    nonzero_prizes_input_file.to_csv(nonzero_output_path, sep='\t', index=False, header=True)
    print(f"Prize file saved for cell line '{cell_line_name}' at: {nonzero_output_path}")

    return gene_to_uniprot

def generate_gold_standard(cell_line_name, model_id, CRISPR_dependency, gene_to_uniprot, threshold: float):
    '''Generate gold standard file for the cell line based on CRISPR dependency and gene to Uniprot mapping.'''
    print(f"Generating gold standard file for cell line '{cell_line_name}'...")
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
    gold_standard_output_path = os.path.join("..", "processed", f"{cell_line_name}_gold_standard.txt")
    gold_standard.to_csv(gold_standard_output_path, sep='\t', index=False, header=False)
    print(f"Gold standard file saved for cell line '{cell_line_name}' at: {gold_standard_output_path}")
    print(f"Number of genes in gold standard: {len(gold_standard)}")

def main():
    print(f"processing cell line: {cell_line_name}")


    # load raw dataset files
    #damaging mutations matrix
    base_dir = os.path.join("..", "raw")  # relative path to raw directory
    damaging_mutations_df = pd.read_csv(os.path.join(base_dir, "OmicsSomaticMutationsMatrixDamaging.csv"),index_col = 0 )
    #omics profiles
    omics_profiles = pd.read_csv(os.path.join(base_dir, "OmicsProfiles.csv"), index_col=0)
    #omics expresion data
    omics_expression = pd.read_csv(os.path.join(base_dir, "OmicsExpressionProteinCodingGenesTPMLogp1.csv"), index_col=0)
    #omics copy number variation data
    omics_cnv = pd.read_csv(os.path.join(base_dir, "OmicsCNGeneWGS.csv"), index_col=0)
    #estimated gene dependency probability based on CRISPR data
    CRISPR_dependency = pd.read_csv(os.path.join(base_dir, "CRISPRGeneDependency.csv"), index_col=0)

    # load uniprot mapping file
    print("loading UniProt mapping file...")
    uniprot_map = pd.read_csv(os.path.join("..", "processed", "DamagingMutations_idMapping_20250718.tsv" ), sep='\t')

    # check if cell line is present in all datasets and proceed
    is_present, model_id = check_cell_line(cell_line_name, omics_profiles, damaging_mutations_df, CRISPR_dependency, omics_expression, omics_cnv)

    if is_present:
        # generate prize input files
        gene_to_uniprot = generate_prize_files(cell_line_name, model_id, damaging_mutations_df, uniprot_map)

        # generate gold standard file
        generate_gold_standard(cell_line_name, model_id, CRISPR_dependency, gene_to_uniprot, dependency_threshold)
        print(f"Processing for cell line '{cell_line_name}' completed successfully.")
    else:
        print(f"Cell line '{cell_line_name}' is not present in all datasets.")
if __name__ == "__main__":
    main()

