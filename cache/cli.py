"""
Downloads the online variants of cache items.

This may be expanded in the future, so only depend on this file as a debugging utility.

For example, `python cache/cli.py a/b.c b.c` would download the file under `a`, `b.c` in `directory`
to the file `b.c`: this is especially useful if headers are provided in `directory.py`
"""

import argparse
from cache.directory import get_cache_item


def parse_args():
    parser = argparse.ArgumentParser(prog="Cache", description="CLI utility for directory.py")
    parser.add_argument("path")
    parser.add_argument("output")

    return parser.parse_args()


def main():
    args = parse_args()
    cache_item = get_cache_item(args.path.split("/"))

    cache_item.download(args.output)


if __name__ == "__main__":
    main()
