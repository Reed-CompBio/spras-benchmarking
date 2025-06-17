from pathlib import Path
import shutil
import pandas

CURRENT_DIR = Path('datasets', 'rn-muscle-skeletal')
PROCESSED_DIR = CURRENT_DIR / 'processed'

def process():
    shutil.copytree(CURRENT_DIR / 'curated', PROCESSED_DIR, dirs_exist_ok=True)

    # TODO: what are the actual last two headers called?
    data = pandas.read_csv(CURRENT_DIR / 'raw' / 'Muscle_Skeletal-Dec2018.tsv',
                           delimiter='\t', header=None,
                           names=["Interactome1", "Interactome2", "Type1",
                                  "Type2", "InteractionType", "Weight",
                                  "Const1", "Const2"])
    data = data.drop(columns=["Type1", "Type2", "InteractionType", "Const1", "Const2"])
    data.insert(3, "Direction", "U")
    data.to_csv(PROCESSED_DIR / 'interactome.tsv', sep='\t', header=False, index=False)

if __name__ == '__main__':
    process()
