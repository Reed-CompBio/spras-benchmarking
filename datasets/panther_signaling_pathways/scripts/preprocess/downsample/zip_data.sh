#!/bin/bash
set -euo pipefail

src="processed"
dst="zip_processed"

find "$src" -name '*.txt' | while read -r f; do
    rel="${f#"$src"/}"                  # <pathway>/interactomes/interactome_X.txt
    pathway=$(echo "$rel" | cut -d/ -f1) # <pathway>  (drops interactomes/)
    base=$(basename "$f" .txt)          # interactome_X
    outdir="$dst/$pathway"
    mkdir -p "$outdir"
    zip -j "$outdir/$base.zip" "$f"
done