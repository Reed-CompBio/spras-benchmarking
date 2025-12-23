import argparse
import gzip
import os
from pathlib import Path
import shutil

from cache.directory import get_cache_item

# https://stackoverflow.com/a/5137509/7589775
dir_path = os.path.dirname(os.path.realpath(__file__))

string_path = Path(dir_path, "string")


def parse_args():
    parser = argparse.ArgumentParser(
        prog="STRING DB Fetcher", description="Downloads specified STRING DB background interactomes from a specific organism."
    )

    parser.add_argument(
        "-i",
        "--id",
        help="""
                        The specified organism ID to use.
                        See https://string-db.org/cgi/download for more info.
                        For example, 9606 is the homo sapiens background interactome.
                        For an example usage, see datasets/diseases's Snakefile.
                        """,
        type=int,
        required=True,
    )

    return parser.parse_args()

def uncompress(source: Path, target: Path):
    """Uncompresses a .gz file"""
    # Uncompressing a .gz file: https://stackoverflow.com/a/44712152/7589775
    with gzip.open(source, "rb") as f_compressed:
        with open(target, "wb") as f_uncompressed:
            shutil.copyfileobj(f_compressed, f_uncompressed)

def main():
    args = parse_args()
    string_path.mkdir(exist_ok=True)

    # We download the links file
    links_file = string_path / f"{args.id}.protein.links.v12.0.txt.gz"
    get_cache_item(["STRING", str(args.id), "links"]).download(links_file)
    uncompress(links_file, links_file.with_suffix("")) # an extra call of with_suffix strips the `.gz` prefix

    # and its associated aliases
    aliases_file = string_path / f"{args.id}.protein.aliases.v12.0.txt.gz"
    get_cache_item(["STRING", str(args.id), "aliases"]).download(aliases_file)
    uncompress(aliases_file, aliases_file.with_suffix(""))

if __name__ == "__main__":
    main()
