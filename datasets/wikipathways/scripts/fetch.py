from pathlib import Path
import os
from cache.directory import get_cache_item
import wikinetworks
import networkx as nx

dir_path = os.path.dirname(os.path.realpath(__file__))

raw_dir = Path(dir_path, "..", "raw")

def runParsePathway(gpml: str) -> nx.DiGraph:
    featureList = wikinetworks.makeFeatureLists(gpml, isFromFile=True)
    featureDFs = wikinetworks.getFeatureDFs(featureList)
    featureDFs['interactDF'] = wikinetworks.mapEndPoints(featureDFs)
    featureDFs['interactDF'] = wikinetworks.processInteractDF(featureDFs, featureList)
    graph = wikinetworks.makeGraph(featureDFs, featureList)
    return graph

def main():
    raw_dir.mkdir(exist_ok=True)

    # TODO: Pass as CLI argument
    pathway = "WP428"

    out_path = raw_dir / f"{pathway}.gpml"
    print(f"Fetching WikiPathways {pathway} GPML...")
    get_cache_item(["WikiPathways", f"{pathway}.gpml"]).download(out_path)

    graph = runParsePathway(out_path.read_text())
    (raw_dir / ".." / "processed").mkdir(exist_ok=True)
    nx.write_graphml(graph, raw_dir / ".." / "processed" / f"{pathway}.graphml")

if __name__ == '__main__':
    main()
