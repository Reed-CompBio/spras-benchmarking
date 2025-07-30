from Bio.KEGG.KGML.KGML_parser import read
from bioservices import UniProt, KEGG
import pandas as pd
from more_itertools import chunked

pathway = read(open("Raw_Data/ko03250.xml", "r"))

# Read in Kegg pathway data and keep only orthologs
entries_data = []
for entry in pathway.entries.values():
    if entry.type == "ortholog":
        entries_data.append({"name": entry.name})
entries_df = pd.DataFrame(entries_data)

# Some orthologs have multiple ko codes in the same row
# The following two lines move all ko codes to individual rows
orthology_ids = entries_df["name"].str.split(" ").explode()
orthology_ids = orthology_ids.apply(lambda x: x.split(":")[1]).tolist()

# Using bioservices KEGG class to map ortholog(ko) codes to human(hsa) codes
k = KEGG()
ko_hsa_map = k.link("hsa", "+".join(orthology_ids))
ko_hsa_dict = {x.split("\t")[0].split(":")[1]: x.split("\t")[1] for x in ko_hsa_map.split("\n")[:-1]}
ko_hsa_df = pd.DataFrame(ko_hsa_dict.items(), columns=["KEGG_Orthology", "HSA"])

# Kegg .get is limited to 10 entries per call
# The following code chunks the hsa list into sets of 10
# then calls the .get function on each which returns kegg api data in string format
hsa_chunked = list(chunked(ko_hsa_df["HSA"].tolist(), 10))
raw_uniprot = []
for entry in hsa_chunked:
    raw_uniprot.append(k.get("+".join(entry)).split("\n///\n\n"))

# Raw Kegg api data is filtered to obtain hsa and uniprot codes for each protein
# Note: Although bioservices .link and .conv return cleaner outputs, they do not support
# one to many relationships at this time.
# Note: bioservices also supplies a parser method for the kegg api but it is also broken at this time.
processed_uniprot = []
for chunk in raw_uniprot:
    for item in chunk:
        item = item.split("\n")
        processed_uniprot.append([(x.strip().split(" ")[1:], "hsa:" + (item[0].split(" " * 7)[1])) for x in item if "UniProt" in x][0])

# Creates a dictionary where uniprot ids are keys and hsa ids are values
hsa_uniprot_dict = {}
for item in processed_uniprot:
    for entry in item[0]:
        hsa_uniprot_dict.update({"up:" + entry: item[1]})

# Creates a dataframe with uniprot and hsa values then merges with ko-hsa dataframe by hsa
hsa_uniprot_map = pd.DataFrame.from_dict(hsa_uniprot_dict.items())
hsa_uniprot_map.columns = ["Uniprot", "HSA"]
final_df = ko_hsa_df.merge(hsa_uniprot_map, on="HSA")
uniprotIDs = final_df["Uniprot"].apply(lambda x: x.split(":")[1]).tolist()

# Filters the combined dataframe to include only rows where the uniprot code is in swissprot
u = UniProt()
tst = u.mapping(fr="UniProtKB", to="UniProtKB-Swiss-Prot", query=",".join(uniprotIDs))
failed_uniprot = pd.Series(list(set(tst["failedIds"]))).apply(lambda x: "up:" + x)

final_df = final_df[~final_df["Uniprot"].isin(failed_uniprot)]
