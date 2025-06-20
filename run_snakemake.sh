#!/bin/sh
# Snakemake does not support importing and orchestrating multiple Snakefiles
# in a way that respects the origin directory.
# Instead, we provide a shell script.

uv tool run snakemake --cores 1 -d datasets/rn-muscle-skeletal -s datasets/rn-muscle-skeletal/Snakefile
uv tool run snakemake --cores 1 -d datasets/yeast-osmotic-stress -s datasets/yeast-osmotic-stress/Snakefile
uv tool run snakemake --cores 1 -d datasets/synthetic-data -s datasets/synthetic-data/Snakefile
