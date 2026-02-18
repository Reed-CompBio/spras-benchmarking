from pathlib import Path
import pandas

egfr_directory = Path(__file__).parent.resolve() / '..'

def main():
    interactome_df = pandas.read_csv(egfr_directory / 'raw' / '9606.protein.links.txt', sep='\t')
    interactome_df['Direction'] = 'U'
    
    (egfr_directory / 'processed').mkdir(exist_ok=True)
    interactome_df.to_csv(egfr_directory / 'processed' / 'interactome.tsv', index=False, header=False, sep='\t')

if __name__ == "__main__":
    main()
