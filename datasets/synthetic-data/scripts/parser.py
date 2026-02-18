import argparse
from pathlib import Path

scripts_directory = Path(__file__).parent.resolve()

def parser():
    parser = argparse.ArgumentParser(prog="PANTHER pathway parser")

    parser.add_argument("pathway", choices=[file.stem for file in (scripts_directory / ".." / "raw" / "pathway-data").iterdir()])

    return parser
