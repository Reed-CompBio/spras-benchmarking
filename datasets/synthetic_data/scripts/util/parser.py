import argparse
from pathlib import Path

from jsonc_parser.parser import JsoncParser

synthetic_directory = Path(__file__).parent.parent.parent.resolve()


def parser():
    parser = argparse.ArgumentParser(prog="PANTHER pathway parser")

    parser.add_argument("pathway", choices=JsoncParser.parse_file(synthetic_directory / "pathways.jsonc"))

    return parser
