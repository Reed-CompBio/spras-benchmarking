import sys
from pathlib import Path
import pandas as pd

processed_dir = Path("processed")
output_file = Path("processed/interactome-stats.tsv")

rows = []

for pathway_dir in sorted(processed_dir.iterdir()):
    if not pathway_dir.is_dir():
        continue
    for interactome_dir in sorted(pathway_dir.glob("interactomes")):
        if not interactome_dir.is_dir():
            continue
        txt_files = list(interactome_dir.glob("interactome_*.txt"))
        if not txt_files:
            continue

        for file in txt_files:
            threshold = file.stem.split("_", 1)[1]

            df = pd.read_csv(file, sep="\t", header=None, usecols=[3], names=["edge_type"])
            counts = df.groupby("edge_type").size()

            row = {
                "pathway": pathway_dir.name,
                "threshold": threshold,
                "total_edges": len(df),
                "undirected_edges": counts.get("U", 0),
                "directed_edges": counts.get("D", 0),
            }
            rows.append(row)


pd.DataFrame(rows).to_csv(output_file, sep="\t", index=False)
print(f"Wrote {len(rows)} rows to {output_file}", file=sys.stderr)