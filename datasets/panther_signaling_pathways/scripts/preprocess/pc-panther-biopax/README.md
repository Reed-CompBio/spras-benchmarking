# panther_pathways

PathwayCommons provides the multi-GB file `pc-biopax.owl`. We need to extract specific pathways from this file.
PaxTools, instead of streaming this XML file, instead opts to load the entire file into memory. Since this is infesable
in any cheap CI system, we instead opt to make this a separate workflow: it takes `pc-biopax.owl`, along with
all PANTHER pathways (TODO: this can be generalized), and generates a new OWL file that contains all PANTHER pathways.

Then, instead of extracting files from the large OWL file above, we use this smaller OWL file in the `../` dataset
where we then split pathways individually.
