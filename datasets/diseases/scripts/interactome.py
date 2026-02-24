from pathlib import Path
import pandas

diseases_path = Path(__file__).parent.parent.resolve()

def main():
    # See /cache/directory.py for information on how this was grabbed.
    # 9606 is the organism code for homo sapiens and the required background interactome of DISEASES.
    string = pandas.read_csv(diseases_path / "raw" / "9606.protein.links.full.txt", sep=" ")
    string = string[["protein1", "protein2", "combined_score"]]

    # Threshold anything above a confidence score of 900 to trim down the background interactome
    string = string[string["combined_score"] > 900]
    # though we still keep the weight afterwards
    (diseases_path / "processed").mkdir(exist_ok=True)
    string.to_csv(diseases_path / "processed" / "string_interactome.tsv", sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()
