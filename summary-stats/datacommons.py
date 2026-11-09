import argparse
import json
from pathlib import Path

'''
script that takes the multiple files in KEGG SPRAS format and checks what nodes from it
are missing from the single big Recon X file
'''
dataset_directory = Path(__file__).parent.parent.resolve()

def get_nodes(sprasfile):
##records set of nodes found in a given metabolic pathway/network file

    nodes = set()
    with open(sprasfile) as fin:
        next(fin)

        for line in fin:
            cols = line.strip().split("\t")
            node1, node2 = cols[0], cols[1]

            nodes.add(node1)
            nodes.add(node2)

    return nodes

def get_edges(sprasfile):
##builds tuple sets of node pairs in KEGG to find in recon

    edges = set()
    with open(sprasfile) as fin:
        next(fin)

        for line in fin:
            cols = line.strip().split("\t")
            node1, node2 = cols[0], cols[1]
            ##cols 1 and 2 in SPRAS format, would be 1 and 3 in SIF
            edge = tuple(sorted((node1, node2)))
            edges.add(edge)

    return edges

def compare(): 
    kegg_dir = dataset_directory / "datasets" / "kegg" / "scripts" / "intermediate" / "spras-kegg-pc-data"
    recon_file = dataset_directory / "datasets" / "reconx" / "spras-recon" / "reconx.hgnc.spras.txt"
    ##currently running this script using SPRAS files, may end up using SIF format if interaction/edge type analysis focus

    recon_edges = get_edges(recon_file)
    recon_nodes = get_nodes(recon_file)

    kegg_edge_not_found = {}
    kegg_node_not_found = {}

    for infile in kegg_dir.glob("*.txt"):
        
        kegg_edges = get_edges(infile)
        kegg_nodes = get_nodes(infile)

        missing_edge = kegg_edges - recon_edges
        ##KEGG - Recon = remaining KEGG edges not found in Recon
        missing_node = kegg_nodes - recon_nodes

        print(f"{infile.name}: {len(kegg_edges)} total edges, with {len(missing_edge)} interactions verbatim not found in Recon file")
        print(f"{infile.name}: {len(kegg_nodes)} total nodes, with {len(missing_node)} nodes verbatim not found in Recon file")
        
        if missing_edge:
            kegg_edge_not_found[infile.name] = missing_edge

        if missing_node:
            kegg_node_not_found[infile.name] = missing_node

    return

##terminal output seems to have majority of interactions and nodes not found in Recon X

if __name__ == "__main__":
    compare()