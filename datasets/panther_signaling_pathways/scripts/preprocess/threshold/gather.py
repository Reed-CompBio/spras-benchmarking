# import glob
# import pandas

# print("here")
# input_files = glob.glob("processed/*/interactomes/overlap_*.tsv")
# print("gathered")
# output_file = "processed/interactome-pathway-overlap.tsv"

# print(input_files)
# print(len(input_files))

# dfs = [pandas.read_csv(f, sep="\t") for f in input_files]
# print("dfs")
# pandas.concat(dfs, ignore_index=True).to_csv(output_file, sep="\t", index=False)
# print("concatted")
# print("outputted")

import glob
import pandas
from concurrent.futures import ThreadPoolExecutor

input_files = glob.glob("processed/*/interactomes/overlap_*.tsv")
print("gathered")
print(len(input_files))
output_file = "processed/interactome-pathway-overlap.tsv"

def read_file(f):
    return pandas.read_csv(f, sep="\t")

with ThreadPoolExecutor() as executor:
    dfs = list(executor.map(read_file, input_files))
print("dfed")

pandas.concat(dfs, ignore_index=True).sort_values(by=["pathway", "threshold"]).to_csv(output_file, sep="\t", index=False)
print("concatted")
print("outputted")