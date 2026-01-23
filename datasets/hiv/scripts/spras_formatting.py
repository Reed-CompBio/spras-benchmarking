from pathlib import Path
import pandas

hiv_directory = Path(__file__).parent.resolve().parent
processed_directory = hiv_directory / "processed"

def main():
    mapping = pandas.read_csv(hiv_directory / 'intermediate' / 'mapping.tsv', sep='\t')
    dict_map = dict(zip(mapping["UniProtKB"], mapping["UniProtKB-ID"]))

    # See prepare.py for the origins of these files.
    prize_05 = pandas.read_csv(hiv_directory / "intermediate" / "prize_05.tsv", sep='\t')
    prize_060 = pandas.read_csv(hiv_directory / "intermediate" / "prize_060.tsv", sep='\t')

    prize_05["Uniprot"] = prize_05["Uniprot"].apply(lambda x: dict_map.get(x))
    prize_060["Uniprot"] = prize_060["Uniprot"].apply(lambda x: dict_map.get(x))
    # We also filter for proteins whose UniProtKB accession numbers no longer exist
    # (usually for being wrongly predicted).
    prize_05 = prize_05[prize_05['Uniprot'].notnull()]
    prize_060 = prize_060[prize_060['Uniprot'].notnull()]

    # Format with SPRAS column names
    prize_05.columns = ["NODEID", "prize"]
    prize_060.columns = ["NODEID", "prize"]

    prize_05.to_csv(processed_directory / "processed_prizes_05.txt", sep="\t", header=True, index=False)
    prize_060.to_csv(processed_directory / "processed_prizes_060.txt", sep="\t", header=True, index=False)

if __name__ == '__main__':
    main()
