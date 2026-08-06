import pandas as pd
from pathlib import Path

dir_path = Path(__file__).parent.resolve()
diseases_path = Path(dir_path, "..")

gs_dir = diseases_path / "GS_files"
prize_dir = diseases_path / "prize_files"
output_csv = diseases_path / "processed" / "gs_prize_overlap.csv"


def main():
    rows = []

    gs_files = sorted(gs_dir.glob("*_GS.txt"))
    print(f'Found {len(gs_files)} gold standard files.')

    for gs_file in gs_files:
        prefix = gs_file.name[: -len("_GS.txt")]
        prize_file = prize_dir / f"{prefix}_prizes.txt"

        if not prize_file.exists():
            print(f'.  Skipping {prefix}: no matching prize file found.')
            continue

        gs_df = pd.read_csv(gs_file, sep="\t", header=None, names=["NodeID"])
        prize_df = pd.read_csv(prize_file, sep="\t")

        gs_nodes = set(gs_df["NodeID"])
        prize_nodes = set(prize_df["NODEID"])

        if len(gs_nodes) == 0:
            print(f'.  Skipping {prefix}: gold standard file is empty.')
            continue

        overlap_nodes = gs_nodes.intersection(prize_nodes)
        overlap_pct = 100 * len(overlap_nodes) / len(gs_nodes)

        rows.append({
            "disease": prefix,
            "overlap_pct": overlap_pct,
            "gs_node_count": len(gs_nodes),
            "prize_node_count": len(prize_nodes),
            "overlap_node_count": len(overlap_nodes),
        })

    result_df = pd.DataFrame(rows).sort_values("disease")
    result_df.to_csv(output_csv, sep="\t", index=False)
    print(f'Wrote overlap results for {len(result_df)} diseases to {output_csv}.')


if __name__ == "__main__":
    main()