'''
Simple stats of given network data
    - number of nodes and edges
        - number of undirected and directed edges respectively
    - average degree
'''

import csv
from pathlib import Path

dataset_directory = Path(__file__).parent.parent.resolve()

def main():
    kegg_dir = dataset_directory / "datasets" / "kegg" / "scripts" / "intermediate" / "spras-kegg-pc-data"
    recon_dir = dataset_directory / "datasets" / "reconx" / "spras-formatted"
    out_dir = dataset_directory / "summary-stats" / "global-sumstats-output"
    out_dir.mkdir(exist_ok=True, parents=True)

    kegg_nets = read_graph(kegg_dir, out_dir, "kegg_sumstats") # KEGG folder vs. Recon X singular file
    recon_net = read_graph(recon_dir, out_dir, "recon_sumstats")    

    get_degree(kegg_nets)
    get_degree(recon_net)

    return

def read_graph(sprasfile, out_dir, out_name):
    '''
    - reads input spras files and returns # of nodes and edges per file
    - returns # of directed & undirected edges within the total edges
    - (should) write summary stats to a .csv or .txt outfile
    '''
    all_edges = []

    rows = ["Name\tNodes\tEdges\tDirected\tUndirected"]
        ## uncomment for .txt
    # rows = ["Name,Nodes,Edges,Directed,Undirected"]
            ## for .csv
        
    for infile in sprasfile.glob("*.txt"):
        # looks through directory provided in main()
        nodes = set()

        edges = []
        directed = []
        undirected = []

        with open(infile) as fin:

            next(fin) #skips counting header line as an edge

            for line in fin:
                row = line.strip().split()
                nodes.add(row[0])
                nodes.add(row[1])
                edges.append(row)
    
                direction = row[3]
                if direction == "D":
                    directed.append(row)
                elif direction == "U":
                    undirected.append(row)

        print(infile.name,':',len(nodes),'nodes and',len(edges),'edges, with', 
              len(directed), 'directed edges', 'and', len(undirected), 'undirected edges') 
            # infile.name gives file name associated w/ the nodes + edges

        print(infile.name,'average degree :', ((2*(len(undirected))+(len(directed))/len(nodes))))
        # average degree of a mixed graph w/ D & U edges: 2E_U + E_D / N

        all_edges.extend(edges)

    ##write sumstats outputs to .csv, .txt, etc
    ##directory for outfile currently exists, write function for individual file outputs 

        rows.append(f"{infile.stem}\t{len(nodes)}\t{len(edges)}\t{len(directed)}\t{len(undirected)}")
        (out_dir / f"{out_name}.txt").write_text("\n".join(rows) + "\n") 
            ## replace .csv w/ .txt if needed, and change delimiters from commas -> \t
        print(f"wrote summary stats to datafile")

    #NOTE: external file doesn't work for recon data (expected). needs solution later

    return all_edges

def get_degree(edgelist):
    # potential use for degree distribution (?)
    # not used for anything rn, ignore

    D = {}

    for e in edgelist:
        D[e[0]] = D.get(e[0],0) + 1
        D[e[1]] = D.get(e[1],0) + 1

    return

if __name__ == '__main__':
    main()