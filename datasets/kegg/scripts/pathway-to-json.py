'''
Script that creates a .json file using the pathway.txt file from Pathway Commons, meant to be used for
[database]-fetchpathway.py. Should create an organized list of pathway name inputs utilizing only the
pathway names per database converted for, from the pathway.txt that features every single PC database.
'''

import json
from pathlib import Path

kegg_directory = Path(__file__).parent.resolve()

def convert():
    raw = (kegg_directory / "raw" / "pathways.txt").read_text().splitlines()
    kegg_pathways = {}

    for line in raw[1:]:  ##skips header
        line = line.strip()
        parts = line.split("\t") ##turns tabs in txt file into individual lines
#        if len(parts) < 2: ##for lines that don't have <= 2 cols, not really an issue for KEGG specifically though
#            continue
        uri, name = parts[0].strip(), parts[1].strip()
        if uri.startswith("http://identifiers.org/kegg"): ##using URI as key phrase
            kegg_pathways[name] = uri

    out = kegg_directory / "intermediate" / "kegg_pathways.json"
    out.parent.mkdir(exist_ok=True, parents=True)
    out.write_text(json.dumps(kegg_pathways, indent=2)) ##json conversion
    print(f"Wrote {len(kegg_pathways)} KEGG pathways to {out}") ##adjust print statement according to later conversions

if __name__ == "__main__":
    convert()