from pathlib import Path
import pandas as pd
import sys
import requests
import time

dir_path = Path(__file__).parent.resolve()

diseases_path = Path(dir_path, "..")
(diseases_path / "processed").mkdir(exist_ok=True, parents=True)

OXO_URL = "https://www.ebi.ac.uk/spot/oxo/api/search?size=500"
BATCH_SIZE = 20

'''
Processes TIGA trait-gene associations into disease-gene associations through namespace mapping. Relies on `gold_standard.csv` to pull the relevant DOIDs.
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
    gold_standard_file = pd.read_csv(diseases_path / "processed" / "gold_standard.csv", sep=",")
    DOIDs = gold_standard_file["diseaseID"].unique().tolist()
    print(f'.  There are {len(DOIDs)} disease ontology IDs from the gold standard.')

    # Querying EBI's OxO ontology mapper to get DOID-to-TIGA namespaces mapped.
    print(f'Querying https://www.ebi.ac.uk/spot/oxo/api/search?size=500 in batches of {BATCH_SIZE}:')
    namespace_mapping = {} # EFO/MONDO -> DOID
    for i in range(0, len(DOIDs), BATCH_SIZE):
        batch = DOIDs[i:i+BATCH_SIZE]
        print(f".  [{i}:{i+BATCH_SIZE}/{len(DOIDs)}]")
        if BATCH_SIZE < 15:
            print(f' --> Mapping {batch}')
        namespace_mapping.update(oxo_search(batch))

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
    tiga_string_df.to_csv(diseases_path / "processed" / "trait_gene_assoc.csv", index=False)

# defined with help from ChatGPT.
def oxo_search(doids, mapping_targets=["EFO", "MONDO"], distance=1,sleep_seconds = 0.01):
    """Map a batch of DOIDs to EFO/MONDO IDs using OxO."""
    
    payload = {"ids": doids, "mappingTarget": mapping_targets, "distance": distance}

    response = requests.post(
        OXO_URL,
        json=payload,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
        timeout=60,
    )

    if response.status_code != 200:
        print(f"WARNING: OxO failed for {doids}: {response.status_code}")
        print(response.text[:500])
        return []

    data = response.json()
    mapped_dict = {}
    for result in data.get("_embedded", {}).get("searchResults", []):
        query_id = result.get("queryId")
        for mapping in result.get("mappingResponseList", []):
            mapped_id = mapping.get("curie")
            if mapped_id is None:
                continue
            mapped_dict[mapped_id] = query_id

    time.sleep(sleep_seconds)
    return mapped_dict    

if __name__ == "__main__":
    main()
