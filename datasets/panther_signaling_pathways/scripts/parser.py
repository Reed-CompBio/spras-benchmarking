import argparse
from pathlib import Path

from jsonc_parser.parser import JsoncParser
import urllib.parse


# TODO: deduplicate from ../Snakefile
def make_file_safe(input_str: str) -> str:
    return urllib.parse.quote(input_str, safe="")


def parser():
    parser = argparse.ArgumentParser(prog="PANTHER pathway parser")


    parser.add_argument("pathway", choices=list(map(make_file_safe, JsoncParser.parse_file("pathways.jsonc"))))

    return parser
