"""
Makes the final SPRAS-necessary adjustments to file prizes.
"""

import argparse
import pandas as pd
from pathlib import Path

def parser():
    parser = argparse.ArgumentParser(
                    prog='PrizeProcessor',
                    description='Makes the final SPRAS-necessary adjustments to file prizes.')

    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    return parser

def main():
    args = parser().parse_args()

    input_file = args.input
    output_file = args.output
    assert input_file is not None, "The input file must be specified."
    assert output_file is not None, "The output file must be specified."

    assert Path(input_file).exists(), "The input file must be on the path!"

    input_df = pd.read_csv(input_file, delimiter='\t', header=None)
    input_df.columns = ['NODEID', 'prize']
    # Ignore _PSUEDONODE
    input_df = input_df[~input_df['NODEID'].str.endswith("_PSEUDONODE")]
    # Sort by NODEID
    input_df.sort_values(['NODEID'])
    # Truncuate to 9 digits + 1 whole digit
    input_df["prize"] = input_df["prize"].apply(lambda x: '%.9f'%x)

    input_df.to_csv(output_file, sep='\t', index=False, header=False)

if __name__ == '__main__':
    main()
