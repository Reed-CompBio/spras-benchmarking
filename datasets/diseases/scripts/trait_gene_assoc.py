from pathlib import Path
import pandas as pd

dir_path = Path(__file__).parent.resolve()

diseases_path = Path(dir_path, "..")
(diseases_path / "intermediate").mkdir(exist_ok=True, parents=True)

'''
Processes TIGA trait-gene associations into disease-gene associations through namespace mapping.
Relies on `gold_standard.csv` to pull the relevant DOIDs.
'''

def main():
    print('Reading TIGA file...')
    # See ../Snakefile for information on this file's origin
    tiga = pd.read_csv(diseases_path / "raw" / "tiga_gene-trait_stats.tsv", sep="\t")
    tiga = tiga[["ensemblId", "efoId","trait", "n_snp", "n_snpw"]]
    print(f'. There are {len(tiga)} trait-gene snp scores from TIGA.')

    # get TIGA traits (they have an underscore '_')
    tiga_ids = tiga["efoId"].unique().tolist()
    print(f'.  There are {len(tiga_ids)} EFO/MONDO/OBA traits.')

    # Get DOIDs from data/gold_standard.csv
    print('Reading gold standard file...')
    gold_standard_file = pd.read_csv(diseases_path / "intermediate" / "gold_standard.csv", sep=",")
    DOIDs = gold_standard_file["diseaseID"].unique().tolist()
    print(f'.  There are {len(DOIDs)} disease ontology IDs from the gold standard.')

    # Using mapping files from OXO2 to map DOID-to-TIGA namespaces.
    # This previously used an API to OXO, which has been taken down.
    # An API for OXO2 is in progress as of June 30, 2026.
    mapping_df = pd.read_csv(diseases_path / "raw" / "oxo2-mappings-inferred.tsv", sep="\t", comment="#")
    # Keep only the DOIDs of interest
    mapping_df = mapping_df[mapping_df["subject_id"].isin(DOIDs)]
    # Keep only EFO and MONDO mappings
    mapping_df = mapping_df[mapping_df["object_id"].str.startswith(("EFO:", "MONDO:"))]
    # Create reverse mapping: EFO/MONDO -> DOID
    namespace_mapping = dict(zip(mapping_df["object_id"], mapping_df["subject_id"]))
    # TIGA uses underscores ('_') instead of colons (':')
    namespace_mapping = {k.replace(':','_'):v for k,v in namespace_mapping.items()}

    num_mapped_DOIDs = len(set(list(namespace_mapping.values())))
    print(f'.  {num_mapped_DOIDs}/{len(DOIDs)} DOIDs are mapped.')

    num_in_tiga = len([x for x in namespace_mapping.keys() if x in tiga_ids])
    print(f'.  {num_in_tiga}/{num_mapped_DOIDs} of mapped DOIDs are in TIGA.')

    # Replace the EFO/MONDO IDs with DOID IDs and drop the ones that don't map.
    print('Mapping TIGA entries to DOID:')
    tiga["efoId"] = tiga["efoId"].map(namespace_mapping)
    tiga = tiga.dropna(subset=["efoId"])
    print(f'.  There are {len(tiga)} trait-gene snp scores after mapping to DOID and dropping N/As.')

    tiga_ids = tiga["efoId"].unique().tolist()
    print(f'.  There are now {len(tiga_ids)} DROID traits in the TIGA dataset. If this is fewer than above, then multiple TIGA ids map to the same DOID.')

    # Mapping ENSG IDs to ENSP IDs through the STRING aliases file
    # given our ENSG and ENSP (non one-to-one!) mapping `string_aliases`.
    string_aliases = pd.read_csv(diseases_path / "raw" / "9606.protein.aliases.txt", sep="\t", usecols=["#string_protein_id", "alias"])
    string_aliases.columns = ["str_id", "ENSP"]
    string_aliases = string_aliases.drop_duplicates()

    # We can create our TIGA-mapped file through merging on the ENSPs
    print('Mapping STRING IDs:')
    tiga_string_df = tiga.merge(string_aliases, left_on="ensemblId", right_on="ENSP", how="inner")
    print(f'.  There are {len(tiga_string_df)} DOID-gene snp scores.')

    ## Write file to inputs.csv
    tiga_string_df.to_csv(diseases_path / "intermediate" / "trait_gene_assoc.csv", index=False)

if __name__ == "__main__":
    main()
