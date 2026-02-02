import argparse
import random
import pandas
from pathlib import Path
import collections

current_directory = Path(__file__).parent.resolve()

# From SPRAS. TODO: import once SPRAS uses pixi
def convert_undirected_to_directed(df: pandas.DataFrame) -> pandas.DataFrame:
    mask = df['Direction'] == 'U'
    new_df = df[mask].copy(deep=True)
    new_df['Interactor1'], new_df['Interactor2'] = new_df['Interactor2'], new_df['Interactor1']
    new_df['Direction'] = 'D'
    df.loc[mask, 'Direction'] = 'D'
    df = pandas.concat([df, new_df], ignore_index=True)
    return df

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

    weight_dicts = collections.OrderedDict(sorted({int(k * 1000): int(v) for k, v in dict(weight_counts.values).items()}.items()))

    # Using a list then creating the set is faster because of the sets rather than the gets.
    print('Creating item samples...')
    full_list: list[int] = []
    curr_v = 0
    for k, v in weight_dicts.items():
        full_list.extend(map(lambda x: x + curr_v, random.sample(range(1, v), round(0.1 * v))))
        curr_v += v
    full_set = set(full_list)

    print('Reading interactome...')
    interactome_df = pandas.read_csv(current_directory / '..' / 'processed' / 'interactome.tsv', header=None, sep='\t',
                                     names=['Interactor1', 'Interactor2', 'Weight', 'Direction'], usecols=[0, 1])
    print('Sampling interactome...')
    interactome_df = interactome_df.iloc[list(full_set)]

    pathway_df = pandas.read_csv(
        current_directory / '..' / 'processed' / pathway / f'{pathway}_gs_edges.txt', sep='\t',
        names=['Interactor1', 'Interactor2', 'Weight', 'Direction'])
    # We consider an undirected edge to be two directed edges
    pathway_df = convert_undirected_to_directed(pathway_df)
    pathway_df = pathway_df[['Interactor1', 'Interactor2']]

    print(f'Merging {pathway} with interactome...')
    pathway_df = pathway_df.merge(interactome_df, how='inner', on=["Interactor1", "Interactor2"])
    
    pathway_df.to_csv(current_directory / 'pathway.tsv', sep='\t', index=False, header=False)
    # TODO: save interactome to a good dir
    # interactome_df.to_csv(current_directory / 'test.tsv', sep='\t', index=False, header=False)

if __name__ == "__main__":
    main()
