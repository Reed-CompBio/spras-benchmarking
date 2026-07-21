import pandas as pd
from pathlib import Path
from tools.normalize.trim_input_nodes import trim_input_nodes_file
from util import build_hgnc_to_ensp

DEPMAP_DIR = Path(__file__).parent.resolve() / ".."
RAW_DIR = DEPMAP_DIR / "raw"
PREPROCESSED_DIR = DEPMAP_DIR / "preprocessed"
PROCESSED_DIR = DEPMAP_DIR / "processed"


def preprocess_cna(df_cna: pd.DataFrame, df_chosen_cell_lines: pd.DataFrame, hgnc_to_ensp: dict[str, str]) -> None:
    """
    Write per-cell-line CNA scores (genes with score in {-2, 2} -> prize of 0.5).
    """

    # restrict to chosen cell lines
    chosen_ccle_ids = set(df_chosen_cell_lines["ccle_id"])
    df_cna = df_cna[[c for c in df_cna.columns if c in chosen_ccle_ids]]
    # map gene ids HGNC_symbol -> ENSP
    df_cna = df_cna.rename(index=hgnc_to_ensp)
    df_cna = df_cna[df_cna.index.isin(hgnc_to_ensp.values())]

    for cell_line in df_cna.columns:
        scores = df_cna[cell_line]
        filtered = scores[scores.isin([-2, 2])].astype(float)
        filtered[:] = 0.5

        cell_dir = PREPROCESSED_DIR / cell_line
        cell_dir.mkdir(parents=True, exist_ok=True)
        filtered.rename("prize").rename_axis("NODEID").to_csv(cell_dir / "cna.csv", sep="\t", header=True)

    print("Finished preprocessing CNA input node data")


def preprocess_somatic_mutations(df_sm: pd.DataFrame, df_chosen_cell_lines: pd.DataFrame, hgnc_to_ensp: dict[str, str]) -> None:
    """Write per-cell-line damaging somatic mutation scores (>= 1.0)."""
    df_sm = df_sm.drop(columns=["SequencingID", "ModelConditionID", "IsDefaultEntryForMC"])
    df_sm = df_sm[df_sm["IsDefaultEntryForModel"] == "Yes"]

    df_sm = df_sm.merge(
        df_chosen_cell_lines, left_on="ModelID", right_on="depmap_id", how="inner"
    ).drop(columns=["depmap_id", "IsDefaultEntryForModel", "ModelID"])

    # strip "(NCBI)" suffix from gene-name columns
    df_sm.columns = df_sm.columns.str.replace(r"\s*\(.*\)$", "", regex=True)
    df_sm = df_sm.set_index("ccle_id")

    # map HGNC -> ENSP, dropping unmapped columns
    mapped_cols = [c for c in df_sm.columns if c in hgnc_to_ensp]
    df_sm = df_sm[mapped_cols].rename(columns=hgnc_to_ensp)

    for cell_line in df_sm.index:
        scores = df_sm.loc[cell_line]
        filtered = scores[scores >= 1.0]

        cell_dir = PREPROCESSED_DIR / cell_line
        cell_dir.mkdir(parents=True, exist_ok=True)
        filtered.rename("prize").rename_axis("NODEID").to_csv(
            cell_dir / "sm.csv", sep="\t", header=True
        )

    print("Finished preprocessing somatic mutation input node data")


def preprocess_tf_activity(df_tfa: pd.DataFrame, df_chosen_cell_lines: pd.DataFrame, hgnc_to_ensp: dict[str, str]) -> None:
    """Write per-cell-line TF activity scores for TFs whose |z-score| exceeds 2.0 (2 standard).

    Z-scores are computed per TF (row-wise) across all chosen cell lines.
    """
    # remap columns: depmap_id -> ccle_id, keeping only chosen cell lines
    depmap_to_ccle = dict(zip(df_chosen_cell_lines["depmap_id"], df_chosen_cell_lines["ccle_id"]))
    df_tfa = df_tfa[[c for c in df_tfa.columns if c in depmap_to_ccle]].rename(columns=depmap_to_ccle)

    # remap index: HGNC_symbol -> ENSP
    df_tfa = df_tfa.rename(index=hgnc_to_ensp)
    df_tfa = df_tfa[df_tfa.index.isin(hgnc_to_ensp.values())]

    # per-TF z-score: subtract row mean, divide by row std
    df_zscore = df_tfa.sub(df_tfa.mean(axis=1), axis=0).div(df_tfa.std(axis=1), axis=0)

    # for each cell line, keep TFs with |z| >= 2sd; score is absolute raw TFA
    for ccle_id in df_zscore.columns:
        mask = df_zscore[ccle_id].abs() >= 2.0
        scores = df_tfa.loc[mask, ccle_id].abs()

        cell_dir = PREPROCESSED_DIR / ccle_id
        cell_dir.mkdir(parents=True, exist_ok=True)
        scores.rename("prize").rename_axis("NODEID").to_csv(cell_dir / "tf.csv", sep="\t", header=True)

    print("Finished preprocessing transcription factor input node data")

def build_spras_node_table(ccle_id: str, df_interactome: pd.DataFrame) -> None:
    """Combine per-cell-line preprocessed inputs into SPRAS node tables (with and without CNA)."""
    cell_dir = PREPROCESSED_DIR / ccle_id

    sm_path = cell_dir / "sm.csv"
    tf_path = cell_dir / "tf.csv"

    # read inputs
    sm = pd.read_csv(sm_path, sep="\t")
    tf = pd.read_csv(tf_path, sep="\t")

    sm["sources"], sm["targets"] = True, False
    tf["sources"], tf["targets"] = False, True

    # process with CNA
    cna_path = cell_dir / "cna.csv"
    cna = pd.read_csv(cna_path, sep="\t")
    cna["sources"], cna["targets"] = True, False

    frames_with_cna = [sm, tf, cna]
    frames_with_cna = [df for df in frames_with_cna if not df.empty]
    combined_with_cna = pd.concat(frames_with_cna, ignore_index=True)

    # collapse duplicates: keep the max prize and OR the role flags
    agg_with_cna = combined_with_cna.groupby("NODEID", as_index=False).agg(
        prize=("prize", "max"),
        sources=("sources", "any"),
        targets=("targets", "any"),
    )
    agg_with_cna["active"] = True
    agg_with_cna = agg_with_cna[["NODEID", "sources", "targets", "prize", "active"]]

    # trim with interactome
    trimmed_with_cna = trim_input_nodes_file(df_interactome, agg_with_cna)

    # process without CNA
    frames_no_cna = [sm, tf]
    frames_no_cna = [df for df in frames_no_cna if not df.empty]
    combined_no_cna = pd.concat(frames_no_cna, ignore_index=True)

    agg_no_cna = combined_no_cna.groupby("NODEID", as_index=False).agg(
        prize=("prize", "max"),
        sources=("sources", "any"),
        targets=("targets", "any"),
    )
    agg_no_cna["active"] = True
    agg_no_cna = agg_no_cna[["NODEID", "sources", "targets", "prize", "active"]]

    trimmed_no_cna = trim_input_nodes_file(df_interactome, agg_no_cna)

    # Check both versions for empty data
    empty_rows = []
    out_dir = PROCESSED_DIR / ccle_id
    out_dir.mkdir(parents=True, exist_ok=True)

    for trimmed, version_label in [(trimmed_with_cna, "input_nodes_with_cna"), (trimmed_no_cna, "input_nodes_without_cna")]:
        n_sources = trimmed["sources"].sum()
        n_targets = trimmed["targets"].sum()
        if n_sources == 0 or n_targets == 0:
            missing = []
            if n_sources == 0:
                missing.append("sources")
            if n_targets == 0:
                missing.append("targets")
            reason = " and ".join(missing)
            empty_rows.append({"ccle_id": ccle_id, "version": version_label, "empty_data": reason})
            print(f"Skipping {ccle_id} ({version_label}): 0 {reason} after trimming")
        else:
            # Save the non-empty version
            filename = "input_nodes.csv" if version_label == "input_nodes_with_cna" else "input_nodes_no_cna.csv"
            trimmed.to_csv(out_dir / filename, sep="\t", index=False, header=True)

    if empty_rows:
        empty_path = PROCESSED_DIR / "empty_cell_lines.csv"
        pd.DataFrame(empty_rows).to_csv(empty_path, sep="\t", index=False, mode="a", header=not empty_path.exists())


def main():
    PREPROCESSED_DIR.mkdir(exist_ok=True)
    PROCESSED_DIR.mkdir(exist_ok=True)

    # load all raw inputs
    df_idmapping = pd.read_csv(RAW_DIR / "9606.protein.aliases.txt", sep="\t")
    df_chosen_cell_lines = pd.read_csv(PROCESSED_DIR / "shared_cell_lines.txt", sep="\t")
    df_cna = pd.read_csv(RAW_DIR / "data_cna_cbioportal_ccle2019.txt", sep="\t", index_col=0)
    df_sm = pd.read_csv(RAW_DIR / "OmicsSomaticMutationsMatrixDamaging_25Q3.csv", index_col=0)
    df_tfa = pd.read_csv(RAW_DIR / "consensus_tfa_march_6.tsv", sep="\t", index_col=0)
    df_interactome = pd.read_csv(PROCESSED_DIR / "interactome.tsv", sep="\t", header=None, names=["Interactor1", "Interactor2", "Weight", "Direction"])

    hgnc_to_ensp = build_hgnc_to_ensp(df_idmapping)

    preprocess_cna(df_cna, df_chosen_cell_lines, hgnc_to_ensp)
    preprocess_somatic_mutations(df_sm, df_chosen_cell_lines, hgnc_to_ensp)
    preprocess_tf_activity(df_tfa, df_chosen_cell_lines, hgnc_to_ensp)

    for ccle_id in df_chosen_cell_lines["ccle_id"]:
        build_spras_node_table(ccle_id=ccle_id, df_interactome=df_interactome)

    print("Finished processing input node data into SPRAS format")


if __name__ == "__main__":
    main()