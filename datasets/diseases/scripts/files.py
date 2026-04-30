import pandas as pd
from pathlib import Path
import sys

dir_path = Path(__file__).parent.resolve()

diseases_path = Path(dir_path, "..")
(diseases_path / "prize_files").mkdir(exist_ok=True, parents=True)
(diseases_path / "GS_files").mkdir(exist_ok=True, parents=True)

'''
Takes the gold standard disease-gene associations from `gold_standard.py` and
the GWAS trait-gene SNP values from `inputs.py` and combines them to generate
SPRAS-formatted prize files and gold standard files.
'''
def main():
    # Gold Standard file from `gold_standard.py`
    GS_string_df = pd.read_csv(diseases_path / "data" / "gold_standard.csv")
    GS_combined_group = GS_string_df.groupby("diseaseName")
    A = set([str(k).lower() for k, v in GS_combined_group])
    print(f'There are {len(A)} diseases from the gold standard.')

    # Inputs file from `inputs.py`
    tiga_string_df = pd.read_csv(diseases_path / "data" / "inputs.csv")
    tiga_group = tiga_string_df.groupby("trait")
    B = set([str(k).lower() for k, v in tiga_group])
    print(f'There are {len(B)} diseases from the TIGA dataset.')

    print(list(A)[:10],list(B)[:10])
    common_groups = A.intersection(B)
    print(len(common_groups),'common groups!!',common_groups)

    # Identify the diseases that are in the gold standard and the inputs.
    GS_string_df = GS_string_df[GS_string_df["diseaseID"].isin(tiga_string_df["id"])]
    GS_combined_group = GS_string_df.groupby("diseaseName")
    GS_combined_dict = {str(k): v for k, v in GS_combined_group}
    print(f'There are {len(GS_combined_dict)} diseases that are in the gold standard and in TIGA.')

    # Filter the SNP dataset for genes in the disease set.

    # UNRESOLVED ISSUE:
    # There is a possibility that the TIGA dataset does not have the GENE_SET_SIZE_MINIMUM
    # requirement for the gold standard diseases because the sources of evidence are different.
    # However, if there are GENE_SET_SIZE_MINIMUM genes in the gold standard (text mining & knowledge based)
    # we should run that disease as a dataset, no matter how many disease genes have SNP scores in TIGA.
    #
    # Resolving this issue is delayed by the bug in gold_standard.py.
    tiga_filtered = tiga_string_df[tiga_string_df["id"].isin(GS_string_df["diseaseID"])]
    tiga_group = tiga_filtered.groupby("trait")
    print(f'There are {len(tiga_group)} TIGA "traits" (diseases)')
    tiga_dict = {k: v for k, v in tiga_group}
    #tiga_count = {x: len(tiga_dict[x]) for x in tiga_dict.keys()}
    #tiga_count_threshold = {k: v for (k, v) in tiga_count.items() if v > 10}

    tiga_threshold = tiga_filtered.loc[tiga_filtered["trait"].isin(list(tiga_dict.keys()))]
    print(f'There are {len(tiga_threshold)} TIGA gene-train pairs')

    tiga_prizes = tiga_threshold.groupby("trait")
    tiga_prize_dict = {str(k): v for k, v in tiga_prizes}

    for disease in tiga_prize_dict.keys():
        df = tiga_prize_dict[disease]
        df = df[["str_id", "n_snpw"]]
        df = df.rename(columns={"str_id": "NODEID", "n_snpw": "prize"})
        df.to_csv(diseases_path / "prize_files" / f"{disease.replace(' ', '_')}_prizes.txt", sep="\t", index=False)

    for disease in GS_combined_dict.keys():
        df = GS_combined_dict[disease]
        df = df[["str_id"]]
        df.to_csv(diseases_path / "GS_files" / f"{disease.replace(' ', '_')}_GS.txt", sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
