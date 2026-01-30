import pandas
from pathlib import Path

current_directory = Path(__file__).parent.resolve()

def main():
    # convert the interactome to SPRAS format
    print('Reading interactome...')
    df = pandas.read_csv(current_directory / '..' / 'raw' / '9606.protein.links.full.v12.0.txt', sep=' ',
                         usecols=['protein1', 'protein2', 'combined_score'])
    df.columns = ['Protein1', 'Protein2', 'Weight']

    # we also want to representatively remove a certain percentage of elements from the interactome,
    # to make sure our interactome downsampling preserves edge weight distributions
    # (we don't care to preserve other major topological properties just yet.)
    # since this file is large, we opt for streaming the interactome for removing edges instead,
    # pre-computing the number of edges that fit into our 'weight' buckets.

    score_max = 1000
    bucket_count = 100
    buckets = list(map(lambda x: x * (score_max / bucket_count), range(bucket_count + 1)))
    print('Counting bucket sizes...')
    df['Weight'].value_counts(bins=buckets, sort=False).to_csv(
        current_directory / '..' / 'processed' / 'bucket-counts.tsv', sep='\t')
    
    print('Processing interactome...')
    df['Weight'] = df['Weight'].div(1000) # scores are from 1-1000: we normalize from 0-1.
    df['Direction'] = 'U'

    print('Saving interactome...')
    df.to_csv(current_directory / '..' / 'processed' / 'interactome.tsv', sep='\t', header=False)



if __name__ == "__main__":
    main()
