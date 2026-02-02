import argparse
import random
import pandas
from pathlib import Path
import collections
from typing import OrderedDict, Optional, NamedTuple
import networkx
import itertools

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

def count_weights() -> OrderedDict[int, int]:
    """Returns an ordered map (lowest to highest weight) from the weight to the number of elements the weight has"""
    weight_counts = pandas.read_csv(current_directory / '..' / 'processed' / 'weight-counts.tsv', sep='\t')
    return collections.OrderedDict(sorted({int(k * 1000): int(v) for k, v in dict(weight_counts.values).items()}.items()))

def sample_interactome(interactome_df: pandas.DataFrame, weight_mapping: OrderedDict[int, int], percentage: float):
    # Using a list then creating the set is faster because of the sets rather than the gets.
    print('Creating item samples...')
    full_list: list[int] = []
    curr_v = 0
    for k, v in weight_mapping.items():
        full_list.extend(map(lambda x: x + curr_v, random.sample(range(1, v), round(percentage * v))))
        curr_v += v
    full_set = set(full_list)

    print('Sampling interactome...')
    return interactome_df.iloc[list(full_set)]

def read_pathway(pathway_name: str) -> pandas.DataFrame:
    """
    Returns the directed-only pathway from a pathway name,
    with columns Interactor1 -> Interactor2.
    """
    pathway_df = pandas.read_csv(
        current_directory / '..' / 'processed' / pathway_name / f'{pathway_name}_gs_edges.txt', sep='\t',
        names=['Interactor1', 'Interactor2', 'Weight', 'Direction'])
    # We consider an undirected edge to be two directed edges
    pathway_df = convert_undirected_to_directed(pathway_df)
    return pathway_df[['Interactor1', 'Interactor2']]

class SourcesTargets(NamedTuple):
    sources: list[str]
    targets: list[str]

def sources_and_targets(pathway_name: str) -> SourcesTargets:
    """
    Returns the sources and targets associated with a particular pathway
    """
    nodes_df = pandas.read_csv(
        current_directory / '..' / 'processed' / pathway_name / f'{pathway_name}_node_prizes.txt', sep='\t',
        usecols=['NODEID', 'sources', 'targets'])
    sources: list[str] = list(nodes_df[nodes_df['sources'] == 'True']['NODEID'])
    targets: list[str] = list(nodes_df[nodes_df['targets'] == 'True']['NODEID'])

    return SourcesTargets(sources, targets)

def main():
    pathway = parser().parse_args().pathway
    # To guarantee similar topologies, we first build up the 10% interactome that preserves the connectedness of our interactome,
    # _then_ we build up the interactome in further 10% increments until we get our desired list of thresholded interactomes.

    print('Reading interactome...')
    interactome_df = pandas.read_csv(current_directory / '..' / 'processed' / 'interactome.tsv', header=None, sep='\t',
                                     names=['Interactor1', 'Interactor2', 'Weight', 'Direction'], usecols=[0, 1])

    # For performance reasons (groupby is quite slow), we sample in the interactome using the pre-computed weight-counts.tsv file
    weight_mapping = count_weights()

    # Get information about the pathway
    pathway_df = read_pathway(pathway)
    sources, targets = sources_and_targets(pathway)

    while attempt_sample(pathway, pathway_df, weight_mapping, interactome_df, sources, targets) is None:
        pass

def attempt_sample(pathway: str,
                   pathway_df: pandas.DataFrame,
                   weight_mapping: OrderedDict[int, int],
                   interactome_df: pandas.DataFrame,
                   sources: list[str], targets: list[str]) -> Optional[tuple[str, str]]:

    interactome_df = sample_interactome(interactome_df, weight_mapping, percentage=0.1)

    print(f'Merging {pathway} with interactome...')
    pathway_df = pathway_df.merge(interactome_df, how='inner', on=["Interactor1", "Interactor2"])
    
    print("Checking for pathway connectedness...")
    has_connection: Optional[tuple[str, str]] = None
    graph = networkx.from_pandas_edgelist(pathway_df, source='Interactor1', target='Interactor2')
    for source, target in itertools.product(sources, targets):
        if graph.has_node(source) and graph.has_node(target) and networkx.has_path(graph, source, target):
            has_connection = (source, target)
            break

    if has_connection:
        print(f'Found connection: {has_connection}')
        pathway_df.to_csv(current_directory / 'pathway.tsv', sep='\t', index=False, header=False)
        # TODO: save interactome to a good dir
        # interactome_df.to_csv(current_directory / 'test.tsv', sep='\t', index=False, header=False)
    return has_connection

if __name__ == "__main__":
    main()
