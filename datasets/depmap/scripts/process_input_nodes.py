# TODO
# process sources
# somatic mutations
# - keep all the genes that are >= 1.0 (they need to be 1.0 and 2.0)
# copy number alteration
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

    # map the ids from HGNC_Symbol to ENSP
    idmapping_df = pd.read_csv(depmap_directory / "raw" / "9606.protein.aliases.txt", sep = "\t")
    df_hgnc = (
        idmapping_df.loc[idmapping_df['source'] == 'Ensembl_HGNC_symbol']
        .drop(columns='source')
        .rename(columns={'#string_protein_id': 'ENSP', 'alias': 'HGNC_symbol'})
        .assign(ENSP=lambda x: x['ENSP'].str.removeprefix('9606.'))
        .reset_index(drop=True)
    )
    print(df_hgnc)


    df_chosen_cell_lines = pd.read_csv(depmap_directory / "raw" / "shared_cell_lines_may_12.txt", sep="\t")
    # print(df_chosen_cell_lines)

    df_copy_number_alteration = pd.read_csv(depmap_directory / "raw" / "data_copy_number_alterations_cbioportal_ccle2019.txt", sep="\t", index_col=0)

    # restrict copy_number_alteration data to the 633 chosen cellines
    # TODO: decide if I want everything to be in ccle ids or in depmap ids; for now do ccle id
    chosen_ccle_ids = set(df_chosen_cell_lines["ccle_id"])
    keep_cols = [c for c in df_copy_number_alteration.columns if c in chosen_ccle_ids]
    df_copy_number_alteration = df_copy_number_alteration[keep_cols]

    # change ids HGNC_symbol from to ENSP for copy number alteration data

    print(df_copy_number_alteration)

    df_copy_number_alteration = (
        df_copy_number_alteration.reset_index()
            .merge(df_hgnc, left_on='Hugo_Symbol', right_on='HGNC_symbol', how='inner')
            .drop(columns=['Hugo_Symbol', 'HGNC_symbol'])
            .set_index('ENSP')
    )

    print(df_copy_number_alteration)


    for i, cell_line in enumerate(df_copy_number_alteration.columns, start=1):
        scores = df_copy_number_alteration[cell_line]
        # Keep rows (genes) with score -2 or 2 (drop 0 and NaN)
        filtered = scores[scores.isin([-2, 2])]
        # TODO: make the actual score .5

        # Replace surviving scores with 0.5
        filtered = filtered.astype(float)
        filtered[:] = 0.5

        cell_dir = depmap_directory / "preprocessed" / cell_line
        cell_dir.mkdir(parents=True, exist_ok=True)

        out_path = cell_dir / "copy_number_alteration.csv"
        filtered.rename("score").to_csv(out_path, index=True, header=True, sep="\t")

    print("Finished preprocessing copy number alteration input node data")

    df_sm = pd.read_csv(depmap_directory / "raw" / "OmicsSomaticMutationsMatrixDamaging_25Q3.csv", sep=",", index_col=0)
    df_sm.drop(columns=["SequencingID", "ModelConditionID", "IsDefaultEntryForMC"], inplace = True)

    # remove the duplicate cell lines; use IsDefaultEntryForModel = Yes
    df_sm = df_sm[df_sm["IsDefaultEntryForModel"] == "Yes"]

    # restrict copy_number_alteration data to the 633 chosen cellines
    df_sm = df_sm.merge(
        df_chosen_cell_lines,
        left_on="ModelID",
        right_on="depmap_id",
        how="inner",
    ).drop(columns=["depmap_id", "IsDefaultEntryForModel"])

    # Move ccle_id next to ModelID
    cols = ["ModelID", "ccle_id"] + [c for c in df_sm.columns if c not in ("ModelID", "ccle_id")]
    df_sm = df_sm[cols]

    # remove (NCBI) from GENE_NAME (NCBI) to get GENE_NAME
    df_sm.columns = df_sm.columns.str.replace(r"\s*\(.*\)$", "", regex=True)

    df_sm = df_sm.set_index("ccle_id").drop(columns="ModelID")

    print(df_sm)

    # id mapping the columns: HGNC to ENSP
    # hgnc_to_ensp = dict(zip(df_hgnc['HGNC_symbol'], df_hgnc['ENSP']))
    # df_sm = df_sm.rename(columns=hgnc_to_ensp)

    hgnc_to_ensp = dict(zip(df_hgnc['HGNC_symbol'], df_hgnc['ENSP']))
    mapped_cols = [c for c in df_sm.columns if c in hgnc_to_ensp]
    df_sm = df_sm[mapped_cols].rename(columns=hgnc_to_ensp)

    print(df_sm)

    for cell_line in df_sm.index:
        scores = df_sm.loc[cell_line]
        filtered = scores[scores >= 1.0]

        cell_dir = depmap_directory / "preprocessed" / cell_line
        cell_dir.mkdir(parents=True, exist_ok=True)

        out_path = cell_dir / "sm.csv"
        filtered.rename("score").rename_axis("gene_name").to_csv(out_path, sep="\t", header=True)

    print("Finished preprocessing somatic mutation input node data")

    df_tfa = pd.read_csv(depmap_directory / "raw" / "consensus_tfa_march_6.tsv", sep="\t", index_col=0)
    df_2sd_tfa = pd.read_csv(depmap_directory / "raw" / "tfs_beyond_2sd_per_cell_line.csv", sep="\t")

    # build a depmap_id -> ccle_id mapping
    depmap_to_ccle = dict(zip(df_chosen_cell_lines["depmap_id"], df_chosen_cell_lines["ccle_id"]))

    # keep only chosen cell lines that actually appear in df_tfa
    keep_cols = [c for c in df_tfa.columns if c in depmap_to_ccle]
    df_tfa = df_tfa[keep_cols].rename(columns=depmap_to_ccle)


    df_2sd_tfa = df_2sd_tfa.merge(
        df_chosen_cell_lines,
        left_on="cell_line",
        right_on="depmap_id",
        how="inner",
    ).drop(columns=["depmap_id"])

    # Move ccle_id next to ModelID
    cols = ["cell_line", "ccle_id"] + [c for c in df_2sd_tfa.columns if c not in ("cell_line", "ccle_id")]
    df_2sd_tfa = df_2sd_tfa[cols]

    print(df_2sd_tfa)

    hgnc_to_ensp = dict(zip(df_hgnc['HGNC_symbol'], df_hgnc['ENSP']))

    def map_tfs(tf_string):
        symbols = [s.strip() for s in tf_string.split(',')]
        mapped = [hgnc_to_ensp[s] for s in symbols if s in hgnc_to_ensp]
        return ', '.join(mapped)

    df_2sd_tfa['tfs'] = df_2sd_tfa['tfs'].apply(map_tfs)
    df_2sd_tfa['n_tfs'] = df_2sd_tfa['tfs'].apply(lambda s: 0 if s == '' else len(s.split(',')))

    print(df_2sd_tfa)

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

    print("Finished preprocessing transcription factor input node data")

    # format the data into SPRAS formatted data
    # all the data
    for _, row in df_chosen_cell_lines.iterrows():
        ccle_id = row["ccle_id"]

        cell_dir = depmap_directory / "preprocessed" / ccle_id

        sm_path = cell_dir / "sm.csv"
        copy_number_alteration_path = cell_dir / "copy_number_alteration.csv"
        tf_path = cell_dir / "tf.csv"

        if not (sm_path.exists() and copy_number_alteration_path.exists() and tf_path.exists()):
            print(f"skipping {ccle_id}: missing one or more input files")
            continue

        # read inputs; assumes each has columns
        sm = pd.read_csv(sm_path, sep="\t")
        copy_number_alteration = pd.read_csv(copy_number_alteration_path, sep="\t")
        print(copy_number_alteration)
        tf = pd.read_csv(tf_path, sep="\t")

        # change the node-id column name across files
        sm = pd.read_csv(sm_path, sep="\t").rename(columns={"gene_name": "NODEID", "score": "prize"})
        copy_number_alteration = pd.read_csv(copy_number_alteration_path, sep="\t").rename(columns={"ENSP": "NODEID", "score": "prize"})
        tf = pd.read_csv(tf_path, sep="\t").rename(columns={"tf": "NODEID", "tfa": "prize"})

        # tag roles
        sm["sources"] = True
        sm["targets"] = False
        copy_number_alteration["sources"] = True
        copy_number_alteration["targets"] = False
        tf["sources"] = False
        tf["targets"] = True

        # stack and collapse duplicates: keep max prize, OR the role flags
        combined = pd.concat([sm, copy_number_alteration, tf], ignore_index=True)
        print(combined)

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

    # TODO: need to also make a version of the input nodes that doesn't include the copy number alteraction
    # without copy_number_alteration

    for _, row in df_chosen_cell_lines.iterrows():
        ccle_id = row["ccle_id"]

        cell_dir = depmap_directory / "preprocessed" / ccle_id

        sm_path = cell_dir / "sm.csv"
        tf_path = cell_dir / "tf.csv"

        if not (sm_path.exists() and copy_number_alteration_path.exists() and tf_path.exists()):
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
        out_path = cell_dir / "input_nodes_no_copy_number_alteration.tsv"
        agg.to_csv(out_path, sep="\t", index=False)

    print("Finished processing input node data into SPRAS format")

if __name__ == "__main__":
    main()