# convert to SPRAS format for both files?

import pandas as pd
from pathlib import Path
from datasets.panther_signaling_pathways.scripts.parser import parser

panther_directory = Path(__file__).parent.parent.resolve()
processed_directory = panther_directory / "processed"
intermediate_directory = panther_directory / "intermediate"
raw_directory = panther_directory / "raw"

aliases = pd.read_csv(raw_directory / "9606.protein.aliases.txt", sep="\t")

ensembl_aliases = aliases[aliases["source"] == "Ensembl_gene"]
ensg_to_ensp = (
    ensembl_aliases.assign(
        ensp=ensembl_aliases["#string_protein_id"].str.removeprefix("9606.")
    )
    .groupby("alias")["ensp"]
    .apply(list)
    .to_dict()
)

uniprot_id_aliases = aliases[aliases["source"] == "UniProt_ID"]
uniprot_id_to_ensp = (
    uniprot_id_aliases.assign(
        uniprot_id=uniprot_id_aliases["#string_protein_id"].str.removeprefix("9606.")
    )
    .groupby("alias")["uniprot_id"]
    .apply(list)
    .to_dict()
)


def process_sources(sources: pd.DataFrame, pathway_file: pd.DataFrame, interactome_file: pd.DataFrame) -> pd.DataFrame:

    sources = sources[["UniProt name"]].copy()

    sources["ENSP"] = sources["UniProt name"].map(lambda e: ",".join(uniprot_id_to_ensp.get(e, [])))
    sources = sources.drop(columns=["UniProt name"])

    # Trim targets to what is in the pathway
    pathway = pd.read_csv(pathway_file, sep="\t")
    pathway_ensps = set(pathway["Node1"]) | set(pathway["Node2"])
    sources = sources[sources["ENSP"].isin(pathway_ensps)].reset_index(drop=True)

    # TODO: trim what is in the interactome(s)
    interactome = pd.read_csv(interactome_file, sep="\t", header=None, names=["Node1", "Node2", "Weight", "Direction"])
    interactome_ensps = set(interactome["Node1"]) | set(interactome["Node2"])
    sources = sources[sources["ENSP"].isin(interactome_ensps)].reset_index(drop=True)

    return sources


def strip_ensp_version(ensp: str) -> str:
    # ENSP00000262053.3 -> ENSP00000262053
    return ensp.split(".")[0]

def process_targets(targets: pd.DataFrame, pathway_file: pd.DataFrame, interactome_file: pd.DataFrame) -> pd.DataFrame:
    targets = targets[["Symbol", "Ensembl", "Protein"]].copy()

    # using the ENSG - ENSP mapping that is in the STRING interactome instead of what is provided in the data
    targets["ENSP"] = targets["Ensembl"].map(lambda e: ",".join(ensg_to_ensp.get(e, [])))
    targets = targets.drop(columns=["Symbol", "Ensembl", "Protein"])

    # Trim targets to what is in the pathway
    pathway = pd.read_csv(pathway_file, sep="\t")
    pathway_ensps = set(pathway["Node1"]) | set(pathway["Node2"])
    targets = targets[targets["ENSP"].isin(pathway_ensps)].reset_index(drop=True)

    # TODO: trim what is in the interactome(s)
    interactome = pd.read_csv(interactome_file, sep="\t", header=None, names=["Node1", "Node2", "Weight", "Direction"])
    interactome_ensps = set(interactome["Node1"]) | set(interactome["Node2"])
    targets = targets[targets["ENSP"].isin(interactome_ensps)].reset_index(drop=True)

    return targets


def spras_format(sources: pd.DataFrame, targets: pd.DataFrame) -> pd.DataFrame:
    source_ensps = set(sources["ENSP"])
    target_ensps = set(targets["ENSP"])
    all_ensps = source_ensps | target_ensps

    spras_df = pd.DataFrame({"NODEID": sorted(all_ensps)})
    spras_df["sources"] = spras_df["NODEID"].isin(source_ensps)
    spras_df["targets"] = spras_df["NODEID"].isin(target_ensps)
    spras_df["prize"] = 1.0
    spras_df["active"] = True

    return spras_df

def main():
    processed_directory.mkdir(exist_ok=True)
    intermediate_directory.mkdir(exist_ok=True)
    pathway = parser().parse_args().pathway

    pathway_folder = intermediate_directory / pathway
    pathway_file = pathway_folder / "pathway.csv"
    interactome_file = processed_directory/ "interactome.tsv"

    targets = pd.read_csv(raw_directory / "Homo_sapiens_TF.tsv", sep="\t")
    processed_targets = process_targets(targets, pathway_file, interactome_file)
    processed_targets.to_csv(pathway_folder / "targets.csv", sep="\t", index=False)

    sources = pd.read_excel(raw_directory / "table_S3_surfaceome.xlsx", sheet_name="in silico surfaceome only", skiprows=1)
    processed_sources = process_sources(sources, pathway_file, interactome_file)
    processed_sources.to_csv(pathway_folder / "sources.csv", sep="\t", index=False)

    spras_df = spras_format(processed_sources, processed_targets)
    spras_df.to_csv(pathway_folder / "input_nodes.csv", sep="\t", index=False)


if __name__ == "__main__":
    main()