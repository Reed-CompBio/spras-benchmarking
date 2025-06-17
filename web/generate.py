from pathlib import Path
import os
import shutil

CURRENT_DIR = Path('web')
DATA_DIR = CURRENT_DIR / 'output' / 'data'

def get_files_rec(directory: Path) -> list[Path]:
    files: list[Path] = []
    for path in directory.iterdir():
        if path.is_dir():
            files += get_files_rec(path)
        else:
            files.append(path)
    return files

if __name__ == '__main__':
    DATA_DIR.mkdir(exist_ok=True)
    if Path('output').exists():
        shutil.move('output', DATA_DIR)
    
    file_markup: list[str] = []
    for file in get_files_rec(DATA_DIR):
        file_loc = file.relative_to(CURRENT_DIR / 'output')
        file_markup.append(f'<ul><a href="{file_loc}">{file_loc}</a></ul>')

    index = (CURRENT_DIR / 'index.html').read_text()
    index = index.replace("(data)", '\n'.join(file_markup))

    (CURRENT_DIR / 'output' / 'index.html').write_text(index)