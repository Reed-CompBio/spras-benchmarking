'''
    Using this script and sif-to-spras.py should yield 84 files in respective output folders, even though
    there should be 87 KEGG pathways present after fetching pathways. The missing files are blank, likely
    from its contents being found in the blacklist. This script excludes processing anything past a blank 
    new line in order to account for removing the extra columns of info + the extended SIF tidbit at the
    end of the original SIF outputs in sif-kegg-pathway-pc-data folder. 

    To double check later: currently omitted files are:
        1. datasets/kegg/scripts/intermediate/sif-kegg-pathway-pc-data/Glycosaminoglycan-biosynthesis---keratan-sulfate.sif
        2. datasets/kegg/scripts/intermediate/sif-kegg-pathway-pc-data/Other-glycan-degradation.sif
        3. datasets/kegg/scripts/intermediate/sif-kegg-pathway-pc-data/Oxidative-phosphorylation.sif
'''

from pathlib import Path

kegg_directory = Path(__file__).parent.resolve()
removing_cols = {"INTERACTION_DATA_SOURCE", "INTERACTION_PUBMED_ID", "PATHWAY_NAMES", "MEDIATOR_IDS"}

def trim_sifs():
    sif_dir = kegg_directory / "intermediate" / "sif-kegg-pathway-pc-data"
    out_dir = kegg_directory / "intermediate" / "clean-sif-kegg-pc-data"
    out_dir.mkdir(exist_ok=True, parents=True)

    for sif in sif_dir.glob("*.sif"):
        lines = []
        for line in sif.read_text().splitlines():
            if not line.strip():
                break
            lines.append(line.split("\t"))

        if not lines:
            continue

        keep = [i for i, col in enumerate(lines[0]) if col not in removing_cols]
        trimmed = ["\t".join(row[i] for i in keep) for row in lines]

        (out_dir / sif.name).write_text("\n".join(trimmed) + "\n")
        print(f"trimmed {sif.name}")


if __name__ == "__main__":
    trim_sifs()