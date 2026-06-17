from pathlib import Path
import pandas as pd


DEPMAP_DIR = Path(__file__).parent.resolve() / ".."
RAW_DIR = DEPMAP_DIR / "raw"
PREPROCESSED_DIR = DEPMAP_DIR / "preprocessed"
PROCESSED_DIR = DEPMAP_DIR / "processed"

FILES = {
    "tf": "tf.csv",
    "sm": "sm.csv",
    "cna": "cna.csv",
}

OUTPUT_FILE = PROCESSED_DIR / "input_node_stats.tsv"

def load_chosen_cell_lines(processed_dir: Path) -> list[str]:
    """Read the list of chosen cell lines from final_cell_lines.txt.

    Assumes one cell line identifier per line, with possible blank lines
    or surrounding whitespace to strip.
    """
    path = processed_dir / "final_cell_lines.csv"
    df = pd.read_csv(path, sep="\t")
    return df["ccle_id"].tolist()


def count_rows(csv_path: Path) -> int:
    """Return the number of NODEID rows in a tf/sm/cna CSV file.

    Returns 0 if the file does not exist, rather than raising, since a
    missing file most likely means that cell line has zero entries for
    that input type.
    """
    if not csv_path.exists():
        return 0
    df = pd.read_csv(csv_path, sep="\t")
    return len(df)


def count_rows(csv_path: Path) -> int:
    """Return the number of NODEID rows in a tf/sm/cna CSV file.

    Returns 0 if the file does not exist, rather than raising, since a
    missing file means that cell line has zero entries for that input type.
    """
    if not csv_path.exists():
        return 0
    df = pd.read_csv(csv_path, sep="\t", header=0)
    return len(df)


def compute_per_cell_line_counts(
    ccle_ids: list[str], preprocessed_dir: Path
) -> pd.DataFrame:
    """Build a DataFrame with one row per cell line and one column per
    input type (tf, sm, cna), containing the row count for each.
    """
    records = []
    for ccle_id in ccle_ids:
        cell_dir = preprocessed_dir / ccle_id
        record = {"ccle_id": ccle_id}
        for input_type, filename in FILES.items():
            record[input_type] = count_rows(cell_dir / filename)
        records.append(record)
    return pd.DataFrame.from_records(records).set_index("ccle_id")


def summarize_min_max(df_counts: pd.DataFrame) -> pd.DataFrame:
    """Given per-cell-line counts, compute min/max for each input type,
    along with which cell line(s) achieve each min/max.
    """
    summary_records = []
    for input_type in FILES:
        col = df_counts[input_type]
        min_val = col.min()
        max_val = col.max()
        # min_cell_lines = ", ".join(col[col == min_val].index)
        # max_cell_lines = ", ".join(col[col == max_val].index)
        summary_records.append(
            {
                "input_type": input_type,
                "min_count": min_val,
                # "min_cell_lines": min_cell_lines,
                "max_count": max_val,
                # "max_cell_lines": max_cell_lines,
            }
        )
    return pd.DataFrame.from_records(summary_records).set_index("input_type")


def write_combined_output(
    df_counts: pd.DataFrame, df_summary: pd.DataFrame, output_path: Path
) -> None:
    """Write per-cell-line counts and the min/max summary into a single
    output file, separated by a header line for readability.
    """
    with open(output_path, "w") as f:
        f.write("# Per-cell-line input node counts\n")
        df_counts.to_csv(f, sep="\t")
        f.write("\n# Min/max summary across cell lines\n")
        df_summary.to_csv(f, sep="\t")


def main() -> None:
    ccle_ids = load_chosen_cell_lines(PROCESSED_DIR)
    df_counts = compute_per_cell_line_counts(ccle_ids, PREPROCESSED_DIR)
    df_summary = summarize_min_max(df_counts)
    write_combined_output(df_counts, df_summary, OUTPUT_FILE)

    print(f"Processed {len(ccle_ids)} cell lines")
    print(df_summary)
    print(f"Full results written to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
