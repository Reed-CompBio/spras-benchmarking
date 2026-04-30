from enum import Enum
from pathlib import Path
import gzip
import shutil
import zipfile

__all__ = ["PostProcessAction"]


def uncompress_gz(source: Path, target: Path):
    """Uncompresses a .gz file"""
    # https://stackoverflow.com/a/44712152/7589775
    with gzip.open(source, "rb") as f_compressed:
        with open(target, "wb") as f_uncompressed:
            shutil.copyfileobj(f_compressed, f_uncompressed)


def uncompress_zip(source: Path, target: Path):
    """Uncompress a .zip file"""
    # https://stackoverflow.com/a/3451150/7589775
    with zipfile.ZipFile(source, "r") as zip_ref:
        zip_ref.extractall(target)


class PostProcessAction(Enum):
    UNCOMPRESS_GZ = "uncompress_gz"
    """Uncompresses a .gz file. The `out_path` here is a file."""
    UNCOMPRESS_ZIP = "uncompress_zip"
    """Uncompresses a .zip file. The `out_path` here is a directory."""

    # NOTE: this can be generalized to multiple outputs and multiple steps, but I'm explicitly avoiding generalization
    # here, as else you run into a general pipeline design pattern when Snakemake can just handle this for us.
    def run_action(self, in_path: Path, out_path: Path):
        match self:
            case self.UNCOMPRESS_GZ:
                uncompress_gz(in_path, out_path)
            case self.UNCOMPRESS_ZIP:
                uncompress_zip(in_path, out_path)
