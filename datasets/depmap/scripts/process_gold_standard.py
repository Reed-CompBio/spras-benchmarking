import pandas as pd
from pathlib import Path
from tools.normalize.trim_list import trim_node_list
from tools.normalize.interactome import get_interactome_nodes
from util import build_hgnc_to_ensp

DEPMAP_DIR = Path(__file__).parent.resolve() / ".."
RAW_DIR = DEPMAP_DIR / "raw"
PREPROCESSED_DIR = DEPMAP_DIR / "preprocessed"
PROCESSED_DIR = DEPMAP_DIR / "processed"

DEPENDENCY_CUTOFF = 0.5

def main():
    PROCESSED_DIR.mkdir(exist_ok=True)

    df_idmapping = pd.read_csv(RAW_DIR / "9606.protein.aliases.txt", sep="\t")

    df_interactome = pd.read_csv(PROCESSED_DIR / "interactome.tsv", sep="\t", header=None, names=["Interactor1", "Interactor2", "Weight", "Direction"])

    df_chosen_cell_lines = pd.read_csv(PROCESSED_DIR / "shared_cell_lines.txt", sep="\t")

    # drop cell line if empty
    df_empty = pd.read_csv(PROCESSED_DIR / "empty_data.csv", sep="\t")
    drop_ccle = set(df_empty["ccle_id"])
    df_chosen_cell_lines = df_chosen_cell_lines[~df_chosen_cell_lines["ccle_id"].isin(drop_ccle)]

    df_gs = pd.read_csv(RAW_DIR / "CRISPRGeneDependency_25Q3.csv", sep=",", index_col=0, header=0)
    df_gs.index.name = "depmap_id"

    df_gs = df_gs.merge(
        df_chosen_cell_lines, left_on="depmap_id", right_on="depmap_id", how="inner"
    ).drop(columns=["depmap_id"])

    # strip "(NCBI)" suffix from gene-name columns
    df_gs.columns = df_gs.columns.str.replace(r"\s*\(.*\)$", "", regex=True)
    df_gs = df_gs.set_index("ccle_id")

    # change the hgnc -> ensp
    hgnc_to_ensp = build_hgnc_to_ensp(df_idmapping)
    mapped_cols = [c for c in df_gs.columns if c in hgnc_to_ensp]
    df_gs = df_gs[mapped_cols].rename(columns=hgnc_to_ensp)

    interactome_nodes = get_interactome_nodes(df_interactome)
    empty_gs_ccle_ids = []

    for ccle_id, row in df_gs.iterrows():
        passing = row[row > DEPENDENCY_CUTOFF].index.tolist()

        trimmed = trim_node_list(interactome_nodes, passing)

        if not trimmed:
            print(f"{ccle_id} is an empty gold standard after interactome trim ({len(passing)} passing before trim)")
            empty_gs_ccle_ids.append(ccle_id)
            continue

        out_dir = PROCESSED_DIR / ccle_id
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "gold_standard.csv"

        pd.DataFrame({"NODEID": sorted(trimmed)}).to_csv(out_path, sep="\t", index=False)

    print("Finished processing gold standard data")

    # write out the list of ccle_ids with empty gs after trimming
    df_final_cell_lines = df_chosen_cell_lines[~df_chosen_cell_lines["ccle_id"].isin(empty_gs_ccle_ids)]
    final_path = PROCESSED_DIR / "final_cell_lines.csv"
    df_final_cell_lines.to_csv(final_path, sep="\t", index=False)

    # append empty gold standard cell lines to empty_data.csv
    if empty_gs_ccle_ids:
        df_empty_gs = pd.DataFrame({
            "ccle_id": empty_gs_ccle_ids,
            "version": "gold_standard",
            "empty_data": "no nodes after interactome trim"
        })
        empty_path = PROCESSED_DIR / "empty_data.csv"
        df_empty_gs.to_csv(empty_path, sep="\t", index=False, mode="a")

    print(f"\n{len(empty_gs_ccle_ids)} cell lines dropped due to empty gs; {len(df_final_cell_lines)} cell lines written to {final_path}")



if __name__ == "__main__":
    main()
