# Python File that takes in a JSON file and converts it into
# SPRAS format

import json
import pandas as pd
import networkx as nx
from pybiomart import Server

server = Server(host="http://ensembl.org")

dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl'])

ensembl_df = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'], filters={'chromosome_name':['1','2']})

with open("files/ResponseNetNetwork.json", "r") as net_json:
    pre_nx = json.load(net_json)
    
RNN = nx.DiGraph()
SN = nx.DiGraph()

for edge in pre_nx["elements"]["edges"]:
    edge = edge["data"]
    
    u = edge["source"]
    v = edge["target"]
    flow = float(edge["flow"])
    
    RNN.add_edge(u, v, flow = flow)
    
with open ("files/benchmark_gamma10.txt", "r") as spras_net:
    for line in spras_net:
        
        tokens = line.strip().split()
        if "Interactor" in tokens:
            continue
        u = tokens[0]
        v = tokens[1]
        flow = float(tokens[2])
        
        SN.add_edge(u, v, flow = flow)

print("Graphs created!") 
print(f"Graph from YL-lab has {RNN.number_of_nodes()} nodes and {RNN.number_of_edges()} edges")
print(f"Graph from SPRAS has {SN.number_of_nodes()} nodes and {SN.number_of_edges()} edges \n")

n_intersection, e_intersection = list(), list()

if RNN.nodes() == SN.nodes():
    print("YL-Lab and SPRAS networks have the same node-set")
    n_intersection = list(RNN.nodes())
else:
    print("YL-Lab and SPRAS networks have different node-sets \n")
    
    n_intersection = list(set(RNN.nodes()) & set(SN.nodes()))
    print(f"Creating list of {len(n_intersection)} nodes in both node-sets:")
    print(n_intersection)
    
if RNN.edges() == SN.edges():
    print("YL-Lab and SPRAS networks have the same edge-set")
    e_intersection = list(SN.edges())
else:
    print("YL-Lab and SPRAS networks have different edge-sets \n")
    
    e_intersection = list(set(RNN.edges()) & set(SN.edges()))
    print(f"Creating list of {len(e_intersection)} edges in both edge-sets:")
    print(e_intersection)