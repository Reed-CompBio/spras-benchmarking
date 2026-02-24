import pandas
from pathlib import Path

egfr_directory = Path(__file__).parent.resolve() / '..'

def main():
    prizes = pandas.read_csv(
        egfr_directory / 'raw' / 'egfr-prizes.txt', sep='\t',
        header=None, names=['NODEID', 'prize']
    )
    prizes = prizes.loc[~prizes['NODEID'].str.endswith('_PSEUDONODE')]
    # TODO: prize: 10 is a magic value.
    prizes = pandas.concat(
        [prizes, pandas.DataFrame({
            'NODEID': ['EGF_HUMAN'],
            'prize': [10],
            'dummy': ['True'],
            'source': ['True']
        })],
        ignore_index=True
    )
    prizes['active'] = 'True'

    prizes.to_csv(egfr_directory / 'processed' / 'prizes-uniprot.txt', index=False, sep='\t')

if __name__ == "__main__":
    main()
