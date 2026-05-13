# TODO
# process sources
# somatic mutations
# - keep all the genes that are >= 1.0 (they need to be 1.0 and 2.0)
# can
# - keep all genes that are -2 or 2s
# tfa
# - keeping all genes that are 2sd

# each cellline will have its own folder of input nodes + gold standard

import pandas as pd
from pathlib import Path

depmap_directory = Path(__file__).parent.resolve() / ".."

def main():

    (depmap_directory / "preprocessed").mkdir(exist_ok=True)
    (depmap_directory / "processed").mkdir(exist_ok=True)

    df_chosen_cell_lines = pd.read_csv(depmap_directory / "raw" / "shared_cell_lines_may_12.txt", sep="\t")
    print(df_chosen_cell_lines)

    df_can = pd.read_csv(depmap_directory / "raw" / "data_can_cbioportal_ccle2019.txt", sep="\t", index_col=0)

    # restrict can data to the 633 chosen cellines
    # TODO: decide if I want everything to be in ccle ids or in depmap ids; for now do ccle id
    chosen_ccle_ids = set(df_chosen_cell_lines["ccle_id"])
    keep_cols = [c for c in df_can.columns if c in chosen_ccle_ids]
    df_can = df_can[keep_cols]

    for i, cell_line in enumerate(df_can.columns, start=1):
        scores = df_can[cell_line]
        # Keep rows (genes) with score -2 or 2 (drop 0 and NaN)
        filtered = scores[scores.isin([-2, 2])]
        # TODO: make the actual score .5

        # Replace surviving scores with 0.5
        filtered = filtered.astype(float)
        filtered[:] = 0.5


        cell_dir = depmap_directory / "preprocessed" / cell_line
        cell_dir.mkdir(parents=True, exist_ok=True)

        out_path = cell_dir / "can.csv"
        filtered.rename("score").to_csv(out_path, index=True, header=True, sep="\t")

    df_sm = pd.read_csv(depmap_directory / "raw" / "OmicsSomaticMutationsMatrixDamaging_25Q3.csv", sep=",", index_col=0)
    df_sm.drop(columns=["SequencingID", "ModelConditionID", "IsDefaultEntryForMC"], inplace = True)

    # remove the duplicate cell lines; use IsDefaultEntryForModel = Yes
    df_sm = df_sm[df_sm["IsDefaultEntryForModel"] == "Yes"]

    # restrict can data to the 633 chosen cellines
    df_sm = df_sm.merge(
        df_chosen_cell_lines,
        left_on="ModelID",
        right_on="depmap_id",
        how="inner",
    ).drop(columns=["depmap_id", "IsDefaultEntryForModel"])

    # Move ccle_id next to ModelID
    cols = ["ModelID", "ccle_id"] + [c for c in df_sm.columns if c not in ("ModelID", "ccle_id")]
    df_sm = df_sm[cols]

    # remove (Unknown) from GENE_NAME (Unknown) to get GENE_NAME
    df_sm.columns = df_sm.columns.str.replace(r"\s*\(.*\)$", "", regex=True)


    df_sm = df_sm.set_index("ccle_id").drop(columns="ModelID")
    # print(df_sm)

    for cell_line in df_sm.index:
        scores = df_sm.loc[cell_line]
        filtered = scores[scores >= 1.0]

        cell_dir = depmap_directory / "preprocessed" / cell_line
        cell_dir.mkdir(parents=True, exist_ok=True)

        out_path = cell_dir / "sm.csv"
        filtered.rename("score").rename_axis("gene_name").to_csv(out_path, sep="\t", header=True)


    df_tfa = pd.read_csv(depmap_directory / "raw" / "consensus_tfa_march_6.tsv", sep="\t", index_col=0)
    df_2sd_tfa = pd.read_csv(depmap_directory / "raw" / "tfs_beyond_2sd_per_cell_line.csv", sep="\t")

    # build a depmap_id -> ccle_id mapping
    depmap_to_ccle = dict(zip(df_chosen_cell_lines["depmap_id"], df_chosen_cell_lines["ccle_id"]))

    # keep only chosen cell lines that actually appear in df_tfa
    keep_cols = [c for c in df_tfa.columns if c in depmap_to_ccle]
    df_tfa = df_tfa[keep_cols].rename(columns=depmap_to_ccle)

    # print(df_tfa)

    df_2sd_tfa = df_2sd_tfa.merge(
        df_chosen_cell_lines,
        left_on="cell_line",
        right_on="depmap_id",
        how="inner",
    ).drop(columns=["depmap_id"])

    # Move ccle_id next to ModelID
    cols = ["cell_line", "ccle_id"] + [c for c in df_2sd_tfa.columns if c not in ("cell_line", "ccle_id")]
    df_2sd_tfa = df_2sd_tfa[cols]

    for _, row in df_2sd_tfa.iterrows():
        ccle_id = row["ccle_id"]
        tfs = [tf.strip() for tf in row["tfs"].split(",")]

        # pull tfs and scores for this cell line
        # TFs are assigned their corresponding absolute value of the TFA scores
        scores = df_tfa[ccle_id].reindex(tfs).abs()

        cell_dir = depmap_directory / "preprocessed" / ccle_id
        cell_dir.mkdir(parents=True, exist_ok=True)

        out_path = cell_dir / "tf.csv"
        scores.rename("tfa").rename_axis("tf").to_csv(out_path, sep="\t", header=True)

    # format the data into SPRAS formatted data

    # all the data
    for _, row in df_chosen_cell_lines.iterrows():
        ccle_id = row["ccle_id"]

        cell_dir = depmap_directory / "preprocessed" / ccle_id

        sm_path = cell_dir / "sm.csv"
        can_path = cell_dir / "can.csv"
        tf_path = cell_dir / "tf.csv"

        if not (sm_path.exists() and can_path.exists() and tf_path.exists()):
            print(f"skipping {ccle_id}: missing one or more input files")
            continue

        # read inputs; assumes each has columns
        sm = pd.read_csv(sm_path, sep="\t")
        can = pd.read_csv(can_path, sep="\t")
        tf = pd.read_csv(tf_path, sep="\t")

        # change the node-id column name across files
        sm = pd.read_csv(sm_path, sep="\t").rename(columns={"gene_name": "NODEID", "score": "prize"})
        can = pd.read_csv(can_path, sep="\t").rename(columns={"Hugo_Symbol": "NODEID", "score": "prize"})
        tf = pd.read_csv(tf_path, sep="\t").rename(columns={"tf": "NODEID", "tfa": "prize"})

        # tag roles
        sm["sources"] = True
        sm["targets"] = False
        can["sources"] = True
        can["targets"] = False
        tf["sources"] = False
        tf["targets"] = True

        # stack and collapse duplicates: keep max prize, OR the role flags
        combined = pd.concat([sm, can, tf], ignore_index=True)

        agg = combined.groupby("NODEID", as_index=False).agg(
            prize=("prize", "max"),
            sources=("sources", "any"),
            targets=("targets", "any"),
        )
        agg["active"] = True

        # reorder to match SPRAS node-table convention
        agg = agg[["NODEID", "prize", "sources", "targets", "active"]]

        cell_dir = depmap_directory / "processed" / ccle_id
        cell_dir.mkdir(parents=True, exist_ok=True)
        out_path = cell_dir / "input_nodes.tsv"
        agg.to_csv(out_path, sep="\t", index=False)

    # TODO: need to also make a version of the input nodes that doesn't include the can
    # without can

    for _, row in df_chosen_cell_lines.iterrows():
        ccle_id = row["ccle_id"]

        cell_dir = depmap_directory / "preprocessed" / ccle_id

        sm_path = cell_dir / "sm.csv"
        tf_path = cell_dir / "tf.csv"

        if not (sm_path.exists() and can_path.exists() and tf_path.exists()):
            print(f"skipping {ccle_id}: missing one or more input files")
            continue

        # read inputs; assumes each has columns
        sm = pd.read_csv(sm_path, sep="\t")
        tf = pd.read_csv(tf_path, sep="\t")

        # change the node-id column name across files
        sm = pd.read_csv(sm_path, sep="\t").rename(columns={"gene_name": "NODEID", "score": "prize"})
        tf = pd.read_csv(tf_path, sep="\t").rename(columns={"tf": "NODEID", "tfa": "prize"})
        # tag roles
        sm["sources"] = True
        sm["targets"] = False
        tf["sources"] = False
        tf["targets"] = True

        # stack and collapse duplicates: keep max prize, OR the role flags
        combined = pd.concat([sm, tf], ignore_index=True)

        # pick the highest prize and if there are duplicates true false, false true for the source target, they are turned into true true with the highest prize
        agg = combined.groupby("NODEID", as_index=False).agg(
            prize=("prize", "max"),
            sources=("sources", "any"),
            targets=("targets", "any"),
        )
        agg["active"] = True

        # reorder to match SPRAS node-table convention
        agg = agg[["NODEID", "prize", "sources", "targets", "active"]]

        cell_dir = depmap_directory / "processed" / ccle_id
        cell_dir.mkdir(parents=True, exist_ok=True)
        out_path = cell_dir / "input_nodes_no_can.tsv"
        agg.to_csv(out_path, sep="\t", index=False)


if __name__ == "__main__":
    main()