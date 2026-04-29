import pandas
from pathlib import Path
from tools.mapping.ensembl_uniprot import idmapping_uniprot_mapping

egfr_directory = Path(__file__).parent.resolve() / ".."


def main():
    # We get specifically the STRING nodes, as the mapping from UniProt overeagerly maps
    string_nodes = pandas.read_csv(
        egfr_directory / "preprocessed" / "ensp" / "interactome.tsv",
        header=None,
        sep="\t",
        names=["Interactor1", "Interactor2", "Weight", "Direction"],
    )
    interactor_series = pandas.concat([string_nodes["Interactor1"], string_nodes["Interactor2"]], ignore_index=True)

    # Re-read the uniprot nodes from `process_gold_standard.py`
    gold_standard_nodes = (egfr_directory / "preprocessed" / "uniprot" / "gold-standard-nodes.txt").read_text().splitlines()
    # and the prizes from `process_prizes.py`
    prizes = pandas.read_csv(egfr_directory / "preprocessed" / "uniprot" / "input-nodes.txt", sep="\t")

    # We grab our UniProt <-> ENSP mapping
    idmapping_df = idmapping_uniprot_mapping(egfr_directory / "raw" / "HUMAN_9606_idmapping_selected.tsv")
    # Trim it with the interactor series
    # Note that the Ensembl_PRO column is ENSP, versus Ensembl which is ENSG.
    idmapping_df = idmapping_df[idmapping_df["Ensembl_PRO"].isin(interactor_series)]

    # Then map the gold standard nodes
    idmapped_gold_standard_nodes_df = pandas.DataFrame(gold_standard_nodes, columns=["UniProtKB-ID"]).merge(
        idmapping_df, on="UniProtKB-ID", how="left"
    )
    idmapped_gold_standard_nodes_df = idmapped_gold_standard_nodes_df.drop(columns=["UniProtKB-ID", "UniProtKB-AC", "Ensembl"])
    idmapped_gold_standard_nodes_df = idmapped_gold_standard_nodes_df[~idmapped_gold_standard_nodes_df["Ensembl_PRO"].isna()]
    gold_standard_nodes = idmapped_gold_standard_nodes_df["Ensembl_PRO"].astype(str).to_list()
    (egfr_directory / "preprocessed" / "ensp" / "gold-standard-nodes.txt").write_text("\n".join(gold_standard_nodes))

    # Then map the input nodes
    idmapped_input_nodes_df = prizes.merge(idmapping_df, left_on="NODEID", right_on="UniProtKB-ID", how="inner")
    idmapped_input_nodes_df = idmapped_input_nodes_df.drop(columns=["UniProtKB-ID", "UniProtKB-AC", "Ensembl", "NODEID"])
    idmapped_input_nodes_df = idmapped_input_nodes_df[~idmapped_input_nodes_df["Ensembl_PRO"].isna()]
    idmapped_input_nodes_df = idmapped_input_nodes_df.rename(columns={"Ensembl_PRO": "NODEID"})
    idmapped_input_nodes_df = idmapped_input_nodes_df[["NODEID", "prize", "active", "dummy", "source"]]
    idmapped_input_nodes_df.to_csv(egfr_directory / "preprocessed" / "ensp" / "input-nodes.txt", sep="\t", index=False)


if __name__ == "__main__":
    main()
