import argparse
import random
import pandas
from pathlib import Path
import collections

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

    weight_dicts = collections.OrderedDict(sorted({int(k): int(v) for k, v in dict(weight_counts.values).items()}.items()))

    full_list = set()
    first_k = 0
    for k, v in weight_dicts.items():
        full_list = full_list.union(set(map(lambda x: x + first_k, random.sample(range(1, v), int(0.1 * v)))))
        first_k += v
 
    interactome_df = pandas.read_csv(current_directory / '..' / 'processed' / 'interactome.tsv', header=None, sep='\t')
    interactome_df = interactome_df.iloc[list(full_list)]
    interactome_df.to_csv(current_directory / 'test.tsv', sep='\t', index=False, header=False)
    

if __name__ == "__main__":
    main()
