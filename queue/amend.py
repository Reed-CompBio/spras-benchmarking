"""
Adds the currently checked out repository commit to the queue,
in the format specified by README.md
"""
import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(prog='Queue Amend',
                                     description='Amends a new commit to the SPRAS-benchmarking queue')
    parser.add_argument('-d', '--dry-run',
                        help='If specified, this will not commit to spras-benchmarking-queue.',
                        action='store_true')

    return parser.parse_args()

if __name__ == '__main__':
    dry_run = parse_args().dry_run

    # https://stackoverflow.com/a/5137509/7589775
    dir_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(dir_path)
    queue_path = Path(dir_path, 'spras-benchmarking-queue')

    # we don't want this to accidentally mess up, so we remove any leftover queue: this should be a small repository regardless
    if queue_path.exists():
        shutil.rmtree(queue_path)

    # (but we still want it!)
    subprocess.run(["git", "clone", "https://github.com/Reed-CompBio/spras-benchmarking-queue", "--depth", "1"])

    # we don't want computer timezone to affect this output, or runtime
    time = datetime.now(timezone.utc)

    # file is of the format {YYYY_MD}/{timestamp}.txt as motivated in README.md
    yyyy_md = time.strftime('%Y-%m')
    timestamp = round(time.timestamp() * 1000)
    timestamped_file = Path(queue_path, yyyy_md, f"{timestamp}.txt")

    # timestamped_file needs the current spras-benchmarking commit
    # but first, let's make sure this is the canonical spras-benchmarking before we commit to some arbitrary queue
    assert subprocess.check_output(["git", "config", "--get", "remote.origin.url"], encoding='utf-8').strip() == \
        'https://github.com/Reed-CompBio/spras-benchmarking.git'
    spras_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], encoding='utf-8').strip()

    assert not timestamped_file.exists(), f"Timestamped file {timestamped_file} already exists?"
    timestamped_file.parent.mkdir(exist_ok=True)
    # newline to meet UNIX convention
    timestamped_file.write_text(spras_commit + '\n')

    os.chdir(queue_path)

    # add and commit our changes
    subprocess.run(['git', 'add', timestamped_file])
    subprocess.run(['git', 'commit', '-m', f'Queueing {spras_commit}'])

    if not dry_run:
        subprocess.run(['git', 'push'])
