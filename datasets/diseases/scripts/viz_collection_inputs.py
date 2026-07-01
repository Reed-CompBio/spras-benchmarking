# 

"""Visualize the distribution of GS nodes and prize nodes as a scatter plot. 
Written with help from Claude.

Matches files in GS_files/ and prize_files/ by their shared prefix
(e.g. DOID_10652_alzheimers_disease_GS.txt <-> DOID_10652_alzheimers_disease_prizes.txt),
counts lines minus 1 (header) in each, and plots GS count (x) vs prize count (y).
"""
from pathlib import Path
import matplotlib.pyplot as plt
from adjustText import adjust_text

dir_path = Path(__file__).parent.resolve()

diseases_path = Path(dir_path, "..")
(diseases_path / "viz").mkdir(exist_ok=True, parents=True)

GS_DIR = Path(diseases_path / "GS_files")
PRIZE_DIR = Path(diseases_path / "prize_files")
OUT_PATH = Path(diseases_path / "viz" / "collection-input-distribution.png")

# count the number of lines (minus the header) of a file
def count_lines(path: Path) -> int:
    with open(path) as f:
        return sum(1 for _ in f) - 1

def main():
    # read the files from the GS and Prize directories
    gs_files = {f.stem.removesuffix("_GS"): f for f in GS_DIR.glob("*.txt")}
    prize_files = {f.stem.removesuffix("_prizes"): f for f in PRIZE_DIR.glob("*.txt")}
    assert len(gs_files)==len(prize_files)

    common = sorted(set(gs_files) & set(prize_files))
    if not common:
        raise SystemExit("No matching GS/prize file pairs found.")

    missing = (set(gs_files) | set(prize_files)) - set(common)
    if missing:
        print(f"Skipped {len(missing)} unmatched file(s): {sorted(missing)}")

    x = [count_lines(gs_files[prefix]) for prefix in common]
    y = [count_lines(prize_files[prefix]) for prefix in common]
    labels = ['_'.join(prefix.split('_')[2:]) for prefix in common]

    plt.figure(figsize=(10,10))
    plt.scatter(x, y, alpha=0.7)

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Nodes in the Gold Standard (high conf disease-gene pairs)")
    plt.ylabel("Nodes with TIGA Prizes (GWAS scores for genes)")
    plt.title("Gold Standard vs. Prize Nodes Per Dataset")
    plt.tight_layout()
    texts = [plt.text(xi, yi, label, fontsize=7) for xi, yi, label in zip(x, y, labels) if xi>800 or yi>800 or yi<=5]
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))
    plt.savefig(OUT_PATH, dpi=150)
    print(f"Plotted {len(common)} pairs -> {OUT_PATH}")


if __name__ == "__main__":
    main()