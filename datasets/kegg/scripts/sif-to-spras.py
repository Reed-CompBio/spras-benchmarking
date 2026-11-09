# script to turn sif files processed from KEGG -> SPRAS (Nodes 1&2, Rank, Dir.) format

from pathlib import Path

kegg_directory = Path(__file__).parent.resolve()
recon_directory = Path(__file__).parent.parent.parent.resolve()

# PC relationship types listed in https://www.pathwaycommons.org/pc2/formats
undirected = {
    "in-complex-with",
    "interacts-with",
    "neighbor-of",
    "reacts-with"
}

directed = {
    "controls-state-change-of",
    "controls-transport-of",
    "controls-phosphorylation-of",
    "controls-expression-of",
    "catalysis-precedes",
    "consumption-controlled-by",
    "controls-production-of",
    "controls-transport-of-chemical",
    "chemical-affects",
    "used-to-produce"
}

# converting the KEGG SIF folder output from sif-cleaner.py into universal SPRAS format
def sif_to_spras():
    #sif_dir = kegg_directory / "intermediate" / "clean-sif-kegg-pc-data"
    ## for kegg
    sif_dir = recon_directory / "reconx"
    ## for recon
    #out_dir = kegg_directory / "intermediate" / "spras-kegg-pc-data"
    ## for kegg
    out_dir = recon_directory / "reconx" / "spras-formatted"
    ## for recon
    out_dir.mkdir(exist_ok=True, parents=True)

    #print(sif_dir)

    for sif in sif_dir.glob("*.sif"):
        rows = ["Node1\tNode2\tRank\tDirection"]
        for line in sif.read_text().splitlines()[1:]:
            cols = line.split("\t")

            node1, relationship, node2 = cols[0], cols[1], cols[2]

            if relationship in undirected:
                direction = "U"
            elif relationship in directed:
                direction = "D"
            
            rows.append(f"{node1}\t{node2}\t1\t{direction}") 
            # currently hard-coded to have unweighted edges (Rank = 1)

        (out_dir / f"{sif.stem}.txt").write_text("\n".join(rows) + "\n")
        print(f"converted {sif.name} to SPRAS format")

if __name__ == "__main__":
    sif_to_spras()