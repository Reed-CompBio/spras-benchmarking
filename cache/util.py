from pathlib import Path
import gzip
import shutil

def uncompress(source: Path, target: Path):
    """Uncompresses a .gz file"""
    # Uncompressing a .gz file: https://stackoverflow.com/a/44712152/7589775
    with gzip.open(source, "rb") as f_compressed:
        with open(target, "wb") as f_uncompressed:
            shutil.copyfileobj(f_compressed, f_uncompressed)
