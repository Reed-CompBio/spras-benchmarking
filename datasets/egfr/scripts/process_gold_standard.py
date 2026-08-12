from pathlib import Path

egfr_directory = Path(__file__).parent.resolve() / ".."


def main():
    # First, we remove all PSUEDONODES (and any duplicates)
    nodes = (egfr_directory / "raw" / "eight-egfr-reference-all.txt").read_text().splitlines()
    nodes = list(set([node for node in nodes if not node.endswith("_PSEUDONODE")]))

    (egfr_directory / "processed").mkdir(exist_ok=True)
    (egfr_directory / "processed" / "uniprot" / "gold-standard-nodes.txt").write_text("\n".join(nodes))


if __name__ == "__main__":
    main()
