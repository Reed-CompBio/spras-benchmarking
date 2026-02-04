"""
Utility CLI for finding pathway critetia from PathwayCommons based on our desired participant count.
This is meant to be interactive for easily examining the available pathways from PathwayCommons over PANTHER
(and perhaps more later!).

See https://www.pathwaycommons.org/pc2/swagger-ui/index.html#/api-controller-v-2 for the API.
"""

import argparse
import requests

from pydantic import BaseModel

SEARCH_URL = "https://www.pathwaycommons.org/pc2/v2/search"

def parser():
    parser = argparse.ArgumentParser(prog="PathwayCommons PANTHER data explorer")

    return parser

# These schemas were manually examined from the API response, and are thus not exhaustive.
class SearchHit(BaseModel):
    uri: str
    name: str
    biopaxClass: str
    numParticipants: int
    numProcesses: int

class SearchResponse(BaseModel):
    numHits: int
    maxHitsPerPage: int
    searchHit: list[SearchHit]

def request(page: int) -> SearchResponse:
    return SearchResponse.model_validate(requests.post(
        'https://www.pathwaycommons.org/pc2/v2/search',
        headers={
            'accept': 'application/json',
            'Content-Type': 'application/json',
        }, json={
            # Indicates a BioPAX pathway
            'q': 'xrefid:P*',
            'type': 'pathway',
            'organism': [
                '9606',
            ],
            'datasource': [
                'panther',
            ],
            'page': page,
        }
    ).json())

def main():
    args = parser().parse_args()

    # TODO: weirdly constructed loop? could be nicer if we use numHits and maxHitsPerPage
    hits: list[SearchHit] = []
    page = 0
    response = request(page)
    print(f"Paginating {page}...")
    while len(response.searchHit) != 0:
        hits.extend(response.searchHit)
        page += 1
        response = request(page)
        print(f"Paginating {page}...")
    
    for hit in hits:
        print(f"({hit.numParticipants}) {hit.name}")

if __name__ == "__main__":
    main()
