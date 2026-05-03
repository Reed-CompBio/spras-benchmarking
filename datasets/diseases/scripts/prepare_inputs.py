import pandas as pd
from pathlib import Path
import sys

dir_path = Path(__file__).parent.resolve()

diseases_path = Path(dir_path, "..")
(diseases_path / "prize_files").mkdir(exist_ok=True, parents=True)
(diseases_path / "GS_files").mkdir(exist_ok=True, parents=True)

'''
Takes the gold standard disease-gene associations from `gold_standard.py` and
the GWAS trait-gene SNP values from `trait_gene_assoc.py` and combines them to generate
SPRAS-formatted prize files and gold standard files.
'''
def main():
    print('Reading gold standard and trait-gene assoc files:')
    # Gold Standard file from `gold_standard.py`
    gold_standard = pd.read_csv(diseases_path / "processed" / "gold_standard.csv")
    gold_standard_grouped = gold_standard.groupby("diseaseID")
    print(f'.  Gold Standard: {len(gold_standard)} disease-gene assocations for {len(gold_standard_grouped)} diseases.')

    # Inputs file from `trait_gene_assoc.py`
    tiga = pd.read_csv(diseases_path / "processed" / "trait_gene_assoc.csv")
    tiga_grouped = tiga.groupby("efoId")
    print(f'.  Trait-gene associations (TIGA): {len(tiga)} trait-gene associations for {len(tiga_grouped)} diseases.')
    
    # assert that the TIGA groups are a subset of the gold_standard
    # (this was by construction in this refactor; might change later)
    combined_groups = set(tiga_grouped.groups.keys()).intersection(set(gold_standard_grouped.groups.keys()))
    assert combined_groups==tiga_grouped.groups.keys()
    print(f'There are {len(combined_groups)} diseases that are in both the gold standard and in TIGA.')
    
    # now, generate the per-disease inputs.
    for diseaseID in tiga_grouped.groups.keys():
        tiga_df = tiga_grouped.get_group(diseaseID)
        gold_standard_df = gold_standard_grouped.get_group(diseaseID)

        # generate file names as DOID_diseaseName.lower()
        disease_name = gold_standard_df["diseaseName"].iloc[0].lower()
        #print(f' {diseaseID} --> {disease_name}')
        filename_prefix = diseaseID.replace(':','_')+'_'+disease_name.replace(' ','_').replace("'","")
        node_filename = f'{filename_prefix}_prizes.txt'
        gs_filename = f'{filename_prefix}_GS.txt'  

        # write the prize file output
        tiga_df = tiga_df[["str_id", "n_snpw"]]
        tiga_df = tiga_df.rename(columns={"str_id": "NODEID", "n_snpw": "prize"})
        #print(f'writing node prizes to {node_filename}')
        tiga_df.to_csv(diseases_path / "prize_files" / node_filename, sep="\t", index=False)

        # write the gold standard output
        #print(f'writing gold standard output to {gs_filename}')
        gold_standard_df.to_csv(diseases_path / "GS_files" / gs_filename, sep="\t", index=False, header=False)
        #print(f'Wrote to {gs_filename}.')
    print(f'Done writing prize and gold standard files for {len(tiga_grouped.groups.keys())} diseases.')


if __name__ == "__main__":
    main()
