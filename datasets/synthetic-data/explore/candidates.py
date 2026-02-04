"""
Utility CLI for finding pathway critetia from PathwayCommons based on our desired participant count.

See https://www.pathwaycommons.org/pc2/swagger-ui/index.html#/api-controller-v-2 for the API.
"""

import argparse
from urllib import request, parse

SEARCH_URL = "https://www.pathwaycommons.org/pc2/v2/search"

def parser():
    parser = argparse.ArgumentParser(prog="PathwayCommons PANTHER data explorer")

    return parser

def main():
    args = parser().parse_args()

    data = parse.urlencode({
        "q": "name:*signaling*",
        "type": "pathway",
        # only Homo Sapiens
        "organism": ["9606"],
        "datasource": ["KEGG"],
        "page": 0
    }).encode()
    req = request.Request(SEARCH_URL, data=data)
    resp = request.urlopen(req)
    print(resp)

if __name__ == "__main__":
    main()
