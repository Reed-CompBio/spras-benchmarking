from pathlib import Path
import pandas
import os

current_directory = Path(os.path.dirname(os.path.realpath(__file__)))
PROCESSED_DIR = current_directory / 'processed'

def process():
    # TODO: what are the actual last two headers called?
    data = pandas.read_csv(current_directory / 'raw' / 'Muscle_Skeletal-Dec2018.tsv',
                           delimiter='\t', header=None,
                           names=["Interactome1", "Interactome2", "Type1",
                                  "Type2", "InteractionType", "Weight",
                                  "Const1", "Const2"])
    data = data.drop(columns=["Type1", "Type2", "InteractionType", "Const1", "Const2"])
    data.insert(3, "Direction", "U")
    data.to_csv(PROCESSED_DIR / 'interactome.tsv', sep='\t', header=False, index=False)

if __name__ == '__main__':
    process()
