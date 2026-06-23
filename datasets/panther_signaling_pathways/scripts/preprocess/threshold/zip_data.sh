#!/bin/bash
set -euo pipefail

src="processed"
dst="zip_processed" # "zip_processed"

find "$src" -name '*.txt' | while read -r f; do
    rel="${f#"$src"/}"              # <pathway>/interactome_X.txt
    pathway=$(dirname "$rel")       # <pathway>
    base=$(basename "$f" .txt)      # interactome_X
    outdir="$dst/$pathway"
    mkdir -p "$outdir"
    zip -j "$outdir/$base.zip" "$f" # -j junks the path, stores filename only
done