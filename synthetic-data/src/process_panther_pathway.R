library(tidyverse)

# Set the path to the directory containing the PANTHER pathway folders using the 'directory' variable
# Set the name of the pathway to process using the 'pathway' variable
directory <- "pathway-data/"
pathways<- c("Apoptosis_signaling", "B_cell_activation", "Beta3_adrenergic_rec", "Cadherin_signaling", "Hedgehog_signaling", "Insulin_IGF", "Interleukin_signaling", "Notch_signaling", "PDGF_signaling", "Ras", "T_cell_activation", "Toll_signaling", "Wnt_signaling", "p38_MAPK", "Nicotinic_acetylchol", "Fas_signaling", "FGF_signaling", "Interferon_gamma_signaling", "JAK_STAT_signaling", "VEGF_signaling")

# Function for processing PANTHER pathways in Extended SIF format into Node/Edge files for benchmarking inputs
process_panther_pathway <- function(pathway_path, pathway_folder) {
  # Read in PANTHER pathway
  panther_pathway <- read.csv(pathway_path, header = FALSE, fill = TRUE, sep = '\t', stringsAsFactors = FALSE)
  
  # Split dataset at the new header
  split_row_index <- which(panther_pathway$V1 == "PARTICIPANT")
  
  # Process edges (rows before the split)
  edges <- panther_pathway[1:(split_row_index - 1), ]
  colnames(edges) <- edges[1, ]  # First row as header
  edges <- edges[-1, ] %>% 
    select(PARTICIPANT_A, INTERACTION_TYPE, PARTICIPANT_B) %>%
    rename(NODE1 = PARTICIPANT_A, NODE2 = PARTICIPANT_B)
  
  # Filter out chebi entries from edges
  edges_rm <- edges %>% filter(str_detect(NODE1, "^chebi:") | str_detect(NODE2, "^chebi:"))
  edges <- edges %>% filter(!str_detect(NODE1, "^chebi:") & !str_detect(NODE2, "^chebi:"))
  
  # Process nodes (rows after the split)
  nodes <- panther_pathway[(split_row_index):nrow(panther_pathway), ] %>%
    select(V1:V5)
  colnames(nodes) <- nodes[1, ]  # First row as header
  nodes <- nodes[-1, ] %>% 
    select(PARTICIPANT, UNIFICATION_XREF) %>%
    rename(NODE = PARTICIPANT, uniprot = UNIFICATION_XREF)
  nodes$uniprot <- str_replace_all(nodes$uniprot, "uniprot:", "")
  
  # Filter out chebi entries from nodes
  nodes_rm <- nodes %>% filter(str_detect(NODE, "^chebi:"))
  nodes <- nodes %>% filter(!str_detect(NODE, "^chebi:"))
  
  write_tsv(edges, paste0(pathway_folder, "EDGES.txt"))
  write_tsv(edges_rm, paste0(pathway_folder, "DEL-EDGES.txt"))
  write_tsv(nodes, paste0(pathway_folder, "NODES.txt"))
  write_tsv(nodes_rm, paste0(pathway_folder, "DEL-NODES.txt"))
  
  # Function for taking a processed Node dataset from process_panther_pathway and checking for TFs
  process_TFs <- function(nodes, pathway_folder) {
    humanTFs <- read_tsv('human-interactome/Homo_sapiens_TF_Uniprot.txt')
    matches <- inner_join(nodes, humanTFs, by = c("uniprot" = "Uniprot_Accession"))
    targets <- matches %>% select(NODE, uniprot)
    
    write_tsv(targets, paste0(pathway_folder, "TARGETS.txt"))
    
    return(targets)
  }
  
  # Function for taking a processed Node dataset from process_panther_pathway and checking for receptors
  process_receptors <- function(nodes, pathway_folder) {
    receptors <- read_tsv('human-interactome/Homo_sapiens_surfaceome.txt') 
    receptors <- receptors %>% 
      select(`UniProt accession`, `Ensembl gene`, `Membranome Almen main-class`) %>%
      filter(`Membranome Almen main-class` == "Receptors")
    matches <- inner_join(nodes, receptors, by = c("uniprot" = "UniProt accession"))
    sources <- matches %>% select(NODE, uniprot)
    
    write_tsv(sources, paste0(pathway_folder, "SOURCES.txt"))
    
    return(sources)
  }
  
  process_prizes <- function(sources, targets, pathway_folder) {
    prizes100 <- rbind(targets, sources) %>% 
      mutate(prizes = 100, active = "true") %>% 
      select(NODE, uniprot, prizes, active)
    
    write_tsv(prizes100, paste0(pathway_folder, "PRIZES-100.txt"))
  }
  
  sources <- process_receptors(nodes, pathway_folder)
  targets <- process_TFs(nodes, pathway_folder)
  prizes <- process_prizes(sources, targets, pathway_folder)
  
  # Return all four datasets as a list
  return(list(edges = edges, edges_rm = edges_rm, nodes = nodes, nodes_rm = nodes_rm, sources = sources, targets = targets, prizes = prizes))
}

for (pathway in pathways) {
  pathway_folder <- paste0(directory, pathway, "/")
  pathway_path <- paste0(directory, pathway, "/", pathway, ".txt")
  result <- process_panther_pathway(pathway_path, pathway_folder)
}

# Access individual datasets
edges <- result$edges
edges_rm <- result$edges_rm
nodes <- result$nodes
nodes_rm <- result$nodes_rm
sources <- result$sources
targets <- result$targets
prizes <- result$prizes