from pathlib import Path
import pandas

hiv_directory = Path(__file__).parent.resolve().parent

def format(prizes: pandas.DataFrame, uniprot_mapping: dict[str, str]):
    prizes["Uniprot"] = prizes["Uniprot"].apply(lambda x: uniprot_mapping.get(x))

    # We also filter for proteins whose UniProtKB accession numbers no longer exist
    # (usually for being wrongly predicted). Older versions of the UniProtKB mapping can be used
    # to preserve these invalid protein codes.
    prizes = prizes[prizes['Uniprot'].notnull()]

    # Format with SPRAS column names
    prizes.columns = ["NODEID", "prize"]

    return prizes

def main():
    # See name_mapping.py for the origins of mapping.tsv
    mapping = pandas.read_csv(hiv_directory / 'intermediate' / 'mapping.tsv', sep='\t')
    uniprot_mapping = dict(zip(mapping["UniProtKB"], mapping["UniProtKB-ID"]))

    # See prepare.py for the origins of these files.
    prize_05 = format(pandas.read_csv(hiv_directory / "intermediate" / "prize_05.tsv", sep='\t'), uniprot_mapping)
    prize_060 = format(pandas.read_csv(hiv_directory / "intermediate" / "prize_060.tsv", sep='\t'), uniprot_mapping)

    prize_05.to_csv(hiv_directory / "processed" / "processed_prizes_05.txt", sep="\t", header=True, index=False)
    prize_060.to_csv(hiv_directory / "processed" / "processed_prizes_060.txt", sep="\t", header=True, index=False)

if __name__ == '__main__':
    main()
