import argparse
import gzip
import os
from pathlib import Path
import shutil

from spras_benchmarking.cache.directory import get_cache_item

# https://stackoverflow.com/a/5137509/7589775
dir_path = os.path.dirname(os.path.realpath(__file__))

string_path = Path(dir_path, 'string')

def parse_args():
    parser = argparse.ArgumentParser(
        prog='STRING DB Fetcher',
        description='Downloads specified STRING DB background interactomes from a specific organism.')

    parser.add_argument('-i', '--id',
                        help="""
                        The specified organism ID to use.
                        See https://string-db.org/cgi/download for more info.
                        For example, 9606 is the homo sapiens background interactome.
                        For an example usage, see datasets/diseases's Snakefile.
                        """,
                        type=int, required=True)

    return parser.parse_args()

def main():
    args = parse_args()
    string_path.mkdir(exist_ok=True)

    output_file = string_path / f"{args.id}.protein.links.v12.0.txt.gz"
    get_cache_item(['STRING', str(args.id)]).download(output_file)

    output_file_uncompressed = string_path / f"{args.id}.protein.links.v12.0.txt"
    # Uncompressing a .gz file: https://stackoverflow.com/a/44712152/7589775
    with gzip.open(output_file, 'rb') as f_compressed:
        with open(output_file_uncompressed, 'wb') as f_uncompressed:
            shutil.copyfileobj(f_compressed, f_uncompressed)

if __name__ == '__main__':
    main()
