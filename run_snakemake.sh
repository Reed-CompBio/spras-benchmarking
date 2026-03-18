#!/usr/bin/env bash

# Snakemake does not support importing and orchestrating multiple Snakefiles
# in a way that respects the origin directory.
# Instead, we provide a shell script.
# We use some error-catching boilerplate from https://sharats.me/posts/shell-script-best-practices/.

# Strict settings
set -o errexit
set -o nounset

# Forcibly use the CWD
cd "$(dirname "$0")"

main() {
    uv run snakemake --cores 1 -d datasets/diseases -s datasets/diseases/Snakefile
    uv run snakemake --cores 1 -d datasets/egfr -s datasets/egfr/Snakefile
}

main "$@"
