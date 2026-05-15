import pandas as pd
from pathlib import Path

DEPMAP_DIR = Path(__file__).parent.resolve() / ".."
RAW_DIR = DEPMAP_DIR / "raw"
PREPROCESSED_DIR = DEPMAP_DIR / "preprocessed"
PROCESSED_DIR = DEPMAP_DIR / "processed"


def main():
    PROCESSED_DIR.mkdir(exist_ok=True)

    df_cna = pd.read_csv(RAW_DIR / "data_cna_cbioportal_ccle2019.txt", sep="\t", index_col=-0)
    df_cna.index.name = None # removes the hugo_symbol index name
    df_somatic = pd.read_csv(RAW_DIR / "OmicsSomaticMutationsMatrixDamaging_25Q3.csv", index_col=0)
    df_somatic = df_somatic.set_index('ModelID').drop(columns=['SequencingID', 'ModelConditionID','IsDefaultEntryForModel', 'IsDefaultEntryForMC'])
    df_crispr_genes = pd.read_csv(RAW_DIR / "CRISPRGeneDependency_25Q3.csv", sep=",", index_col=0, header=0)
    df_tfs = pd.read_csv(RAW_DIR / "tfs_beyond_2sd_per_cell_line.csv", sep="\t")

    # ccle id -> depmap id; comes from CCLE 2019
    cell_line_annotations = pd.read_csv(RAW_DIR / "Cell_lines_annotations_20181226_ccle2019.txt", sep="\t", header=0)
    cell_line_annotations = cell_line_annotations[["CCLE_ID", "depMapID"]].dropna(subset=["CCLE_ID", "depMapID"])

    # depmap id -> ccle id; comes from DepMap 25Q3
    model_annotations = pd.read_csv(RAW_DIR / "Model_25Q3.csv")
    model_annotations = model_annotations[["ModelID", "CCLEName"]]

    # converting CCLE cell line ids to depmap cell line ids using cell_line_annotations (from ccle 2019)
    cell_line_annotations_ccle_to_depmap_map = dict(zip(cell_line_annotations["CCLE_ID"], cell_line_annotations["depMapID"]))
    df_cna = df_cna.rename(columns=cell_line_annotations_ccle_to_depmap_map)

    # removing columns from df_cna with no matching depMapID
    matched_cols = [col for col in df_cna.columns if col.startswith("ACH-")]
    df_cna = df_cna[matched_cols]

    # remove any cell lines with no copy number alterations (no 2s or -2s)
    cna = df_cna.isin([2, -2]).any(axis=0)
    df_cna = df_cna.loc[:, cna]

    # remove any cell lines with no somatic mutations (all < 1.0)

    has_mutation = (df_somatic >= 1.0).any(axis=1)
    df_somatic = df_somatic.loc[has_mutation, :]

    # remove any cell lines with no tfs
    df_tfs = df_tfs[df_tfs["n_tfs"] >= 1]

    # remove any cell lines with no crispr genes >=0.5
    has_crispr = (df_crispr_genes >= 0.5).any(axis=1)
    df_crispr_genes = df_crispr_genes.loc[has_crispr, :]

    # find all cell lines that match between all of the data
    # all are depmap ids
    cna_cell_lines = set(df_cna.columns)
    somatic_cell_lines = set(df_somatic.index)
    tfs_cell_lines = set(df_tfs["cell_line"])
    crispr_cell_lines = set(df_crispr_genes.index)

    shared_cell_lines = cna_cell_lines & somatic_cell_lines & tfs_cell_lines & crispr_cell_lines

    # Check which cell line annotations and model.csv ids don't having matching depmap ids and ccle ids
    # all of those will be marked as problematic

    # same CCLE ID, different depmap ID
    merged_on_ccle = cell_line_annotations.merge(model_annotations, left_on="CCLE_ID", right_on="CCLEName", how="inner")
    ccle_depmap_mismatch = merged_on_ccle[merged_on_ccle["depMapID"] != merged_on_ccle["ModelID"]]

    # same depmap ID, different CCLE ID
    merged_on_depmap = cell_line_annotations.merge(model_annotations, left_on="depMapID", right_on="ModelID", how="inner")
    depmap_ccle_mismatch = merged_on_depmap[merged_on_depmap["CCLE_ID"] != merged_on_depmap["CCLEName"]]

    # union of all problematic depmap IDs
    problematic_depmap_ids = set(ccle_depmap_mismatch["depMapID"]) | set(ccle_depmap_mismatch["ModelID"]) | set(depmap_ccle_mismatch["depMapID"])

    # remove problematic dep map ids from total number of shared_cell_lines
    shared_cell_lines = shared_cell_lines - problematic_depmap_ids

    shared_df = (
        model_annotations[model_annotations["ModelID"].isin(shared_cell_lines)]
        .sort_values("ModelID")
        .rename(columns={"ModelID": "depmap_id", "CCLEName": "ccle_id"})
    )
    shared_df.to_csv(PROCESSED_DIR / "shared_cell_lines.txt", index=False, sep= "\t")
    print(f"Finished getting shared cell lines: {len(shared_cell_lines)}")


if __name__ == "__main__":
    main()