#!/bin/sh
# Snakemake does not support importing and orchestrating multiple Snakefiles
# in a way that respects the origin directory.
# Instead, we provide a shell script.

snakemake -d datasets/rn-muscle-skeletal -s datasets/rn-muscle-skeletal/Snakefile
snakemake -d datasets/yeast-osmotic-stress -s datasets/yeast-osmotic-stress/Snakefile
