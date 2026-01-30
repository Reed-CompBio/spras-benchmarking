import argparse
import pandas
from pathlib import Path

current_directory = Path(__file__).parent.resolve()

def parser():
    parser = argparse.ArgumentParser(prog='PANTHER pathway parser')

    parser.add_argument('pathway', choices=[file.stem for file in (current_directory / '..' / 'raw' / 'pathway-data').iterdir()])

    return parser

def main():
    pathway = parser().parse_args().pathway

    # To guarantee similar topologies, we first build up the 10% interactome that preserves the connectedness of our interactome,
    # _then_ we build up the interactome in further 10% increments until we get our desired list of thresholded interactomes.

    # For performance reasons (groupby is quite slow), we stream in the interactome using the pre-computed weight-counts.tsv file
    weight_counts = pandas.read_csv(current_directory / '..' / 'processed' / 'weight-counts.tsv', sep='\t')
    total_rows = weight_counts['count'].sum()

    interactome_path = current_directory / '..' / 'processed' / 'interactome.tsv'

    # 

if __name__ == "__main__":
    main()
