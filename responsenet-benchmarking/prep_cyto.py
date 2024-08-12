# Python file that prepares input files for cytoscape analysis

import json
import pandas as pd
import networkx as nx
from pybiomart import Server
import argparse

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
    
    with open(RNN_file, "r") as net_json:
        pre_nx = json.load(net_json)
    
    for edge in pre_nx["elements"]["edges"]:
        edge = edge["data"]
        
        u = get_gene_name(edge["source"])
        v = get_gene_name(edge["target"])
        
        if u not in RNN:
            RNN.add_node(u, ensembl = edge["source"], source_set = False, target_set = False)
        if v not in RNN:
            RNN.add_node(v, ensembl = edge["target"], source_set = False, target_set = False)
        
        flow = float(edge["flow"])
        
        RNN.add_edge(u, v, flow = flow)

    for s_node in nx.all_neighbors(RNN, 'S'):
        RNN.nodes[s_node]["source_set"] = True
    for t_node in nx.all_neighbors(RNN, 'T'):
        RNN.nodes[t_node]["target_set"] = True
        
    with open (SN_file, "r") as spras_net:
        for line in spras_net:
            
            tokens = line.strip().split()
            if "Interactor" in tokens:
                continue
            
            u = get_gene_name(tokens[0])
            v = get_gene_name(tokens[1])
            
            if u not in SN:
                if u == 'source':
                    u = 'S'
                    SN.add_node(u, ensembl = 'S', source_set = False, target_set = False)
                else:
                    SN.add_node(u, ensembl = tokens[0], source_set = False, target_set = False)
                
            if v not in SN:
                if v == 'target':
                    v = 'T'
                    SN.add_node(v, ensembl = 'T', source_set = False, target_set = False)
                else:
                    SN.add_node(v, ensembl = tokens[1], source_set = False, target_set = False)
            
            flow = float(tokens[2])
            
            SN.add_edge(u, v, flow = flow)
            
    for s_node in nx.all_neighbors(SN, 'S'):
        SN.nodes[s_node]["source_set"] = True
    for t_node in nx.all_neighbors(SN, 'T'):
        SN.nodes[t_node]["target_set"] = True
    
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
        
    return n_intersection, e_intersection
        
def write_cytoscape(out_path, json, spras, n_int, e_int):
    RNN_eout = json + "_edges_out" + ".txt"
    RNN_nout = json + "nodes_out" + ".txt"
    with open(RNN_eout, 'w') as RNN_ef:
        RNN_ef.write("Interactor1" + '\t' + "Interactor2" + '\t' + "Flows" + "\t" + "Intersection" + "\n")
        
        for edge in RNN.edges():
            e_intersection = edge in e_int
            RNN_ef.write(edge[0] + "\t" + edge[1] + "\t" + str(RNN[edge[0]][edge[1]]['flow']) + '\t' + str(e_intersection) + "\n")
    
    with open(RNN_nout, 'w') as RNN_nf:
        RNN_nf.write("Node Name" + '\t' + "Set type" + '\t' + "Intersection" +'\n')
        for node in RNN.nodes():
            n_intersection = node in n_int
            if RNN.nodes[node]['source_set']:
                set_type = 1
            elif RNN.nodes[node]['target_set']:
                set_type = 2
            else:
                set_type = 0
            RNN_nf.write(node + '\t' + str(set_type) +'\t'+str(n_intersection)+ '\n')
    
    SN_eout = spras + "_edges_out" + ".txt"
    SN_nout = spras + "nodes_out" + ".txt"
    with open(SN_eout, 'w') as SN_ef:
        SN_ef.write("Interactor1" + '\t' + "Interactor2" + '\t' + "Flows" + "\t" + "Intersection" + "\n")
        
        for edge in SN.edges():
            e_intersection = edge in e_int
            SN_ef.write(edge[0] + "\t" + edge[1] + "\t" + str(SN[edge[0]][edge[1]]['flow']) + '\t' + str(e_intersection) + "\n")
    
    with open(SN_nout, 'w') as SN_nf:
        SN_nf.write("Node Name" + '\t' + "Set type" + '\t' + "Intersection" +'\n')
        for node in SN.nodes():
            n_intersection = node in n_int
            if SN.nodes[node]['source_set']:
                set_type = 1
            elif SN.nodes[node]['target_set']:
                set_type = 2
            else:
                set_type = 0
            SN_nf.write(node + '\t' + str(set_type) +'\t'+str(n_intersection)+ '\n')

        
    
def main(args):
    create_graphs(args.json, args.spras)
    n_int, e_int = compare_graphs()
    write_cytoscape(args.output, args.json, args.spras, n_int, e_int)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--json',
                        help='ResponseNet V.3 Output, file should be in .json format',
                        type=str,
                        required=True)
    parser.add_argument('--spras',
                        help='SPRAS output, file should be in .txt format',
                        type=str,
                        required=True)
    parser.add_argument('--output',
                        help='Filepath to where output file should be written',
                        type=str,
                        required=True)
    
    args = parser.parse_args()
    print(args)

    main(args)