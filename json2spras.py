# Python File that takes in a JSON file and converts it into
# SPRAS format

import json
import pandas as pd
import networkx as nx
from pybiomart import Server

server = Server(host="http://ensembl.org")

dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl'])

ensembl_df = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'], filters={})

RNN = nx.DiGraph()
SN = nx.DiGraph()

def get_gene_name(node):
    translation = ensembl_df[ensembl_df["Gene stable ID"] == node]
    if not translation.empty:
        gene_name = translation.iloc[0]["Gene name"]
    else:
        gene_name = node
        print('not found')
    return gene_name

def create_graphs(RNN_file, SN_file):
    global RNN
    global SN
    
    with open("files/ResponseNetNetwork.json", "r") as net_json:
        pre_nx = json.load(net_json)
    
    for edge in pre_nx["elements"]["edges"]:
        edge = edge["data"]
        
        u = get_gene_name(edge["source"])
        v = get_gene_name(edge["target"])
        
        if u not in RNN:
            RNN.add_node(u, ensembl = edge["source"])
        if v not in RNN:
            RNN.add_node(v, ensembl = edge["target"])
        
        flow = float(edge["flow"])
        
        RNN.add_edge(u, v, flow = flow)
        
    with open ("files/benchmark_gamma10.txt", "r") as spras_net:
        for line in spras_net:
            
            tokens = line.strip().split()
            if "Interactor" in tokens:
                continue
            u = get_gene_name(tokens[0])
            v = get_gene_name(tokens[1])
            
            if u not in SN:
                SN.add_node(u, ensembl = tokens[0])
            if v not in SN:
                SN.add_node(v, ensembl = tokens[1])
            
            flow = float(tokens[2])
            
            SN.add_edge(u, v, flow = flow)
    
    print("Graphs created!") 
    print(f"Graph from YL-lab has {RNN.number_of_nodes()} nodes and {RNN.number_of_edges()} edges")
    print(f"Graph from SPRAS has {SN.number_of_nodes()} nodes and {SN.number_of_edges()} edges \n")
    
def compare_graphs():
    global RNN
    global SN
    
    n_intersection, e_intersection = list(), list()
    # Deleting this soon, cause its doo doo
    if RNN.nodes() == SN.nodes():
        print("YL-Lab and SPRAS networks have the same node-set")
        n_intersection = list(RNN.nodes())
    else:
        print("YL-Lab and SPRAS networks have different node-sets \n")
        
        n_intersection = list(set(RNN.nodes())&set(SN.nodes()))
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
        
        
        
        
    
def main():
    create_graphs("foo", "bar")
    compare_graphs()
    
if __name__ == "__main__":
    main()