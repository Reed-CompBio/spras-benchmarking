import pandas
from pathlib import Path

egfr_directory = Path(__file__).parent.resolve() / '..'

def main():
    prizes = pandas.read_csv(
        egfr_directory / 'raw' / 'egfr-prizes.txt', sep='\t',
        header=None, names=['NODEID', 'prize']
    )
    prizes['active'] = 'True'
    prizes['dummy'] = 'True'

    prizes.to_csv(egfr_directory / 'processed' / 'prizes-uniprot.txt', index=False, sep='\t')

if __name__ == "__main__":
    main()
