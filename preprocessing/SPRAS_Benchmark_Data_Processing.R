library(tidyverse)
# install.packages("httr")
library(httr)
# Set working directory so that pathway_path does not have any parent directories
# Based on GitHub file structure this should be the specific folder for each PANTHER pathway
panther_file <- "Nicotinic_acetylchol"
setwd(paste0('~/Desktop/GitHub_Repos/spras-benchmarking/preprocessing/', panther_file))
pathway_path <- paste0(panther_file, ".txt")

# Function for processing PANTHER pathways in Extended SIF format into Node/Edge files for benchmarking inputs
process_panther_pathway <- function(pathway_path) {
  # Read in PANTHER pathway
  panther_pathway <- read.csv(pathway_path, header = FALSE, fill = TRUE, sep = '\t', stringsAsFactors = FALSE)
  
  # Split dataset at the new header
  split_row_index <- which(panther_pathway$V1 == "PARTICIPANT")
  
  # Process edges (rows before the split)
  edges <- panther_pathway[1:(split_row_index - 1), ]
  colnames(edges) <- edges[1, ]  # First row as header
  edges <- edges[-1, ] %>% 
    select(PARTICIPANT_A, PARTICIPANT_B) %>%
    rename(id1 = PARTICIPANT_A, id2 = PARTICIPANT_B)
  
  # Filter out chebi entries from edges
  edges_rm <- edges %>% filter(str_detect(id1, "^chebi:") | str_detect(id2, "^chebi:"))
  edges <- edges %>% filter(!str_detect(id1, "^chebi:") & !str_detect(id2, "^chebi:"))
  
  # Process nodes (rows after the split)
  nodes <- panther_pathway[(split_row_index):nrow(panther_pathway), ] %>%
    select(V1:V5)
  colnames(nodes) <- nodes[1, ]  # First row as header
  nodes <- nodes[-1, ] %>% 
    select(PARTICIPANT, UNIFICATION_XREF) %>%
    rename(id = PARTICIPANT, uniprot = UNIFICATION_XREF)
  nodes$uniprot <- str_replace_all(nodes$uniprot, "uniprot:", "")
  
  # Filter out chebi entries from nodes
  nodes_rm <- nodes %>% filter(str_detect(id, "^chebi:"))
  nodes <- nodes %>% filter(!str_detect(id, "^chebi:"))
  
  # Check if namespace job is ready
  isJobReady <- function(jobId) {
    pollingInterval = 5
    nTries = 20
    for (i in 1:nTries) {
      url <- paste("https://rest.uniprot.org/idmapping/status/", jobId, sep = "")
      r <- GET(url = url, accept_json())
      # status <- content(r, as = "parsed")
      status <- r[["status_code"]]
      print(paste("Status", status))
      if (status == 200) {
        return(TRUE)
      }
      else {
        return (FALSE)
      }
      Sys.sleep(pollingInterval)
    }
    return(FALSE)
  }
  
  # Get URL of result
  getResultsURL <- function(redirectURL) {
    if (grepl("/idmapping/results/", redirectURL, fixed = TRUE)) {
      url <- gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
    } else {
      url <- gsub("/results/", "/results/stream/", redirectURL)
    }
  }
  
  ### Get STRING IDs from UniProt IDs
  # Set input for namespace mapper
  uniprot_ids <- paste(nodes$uniprot, collapse = ",")
  files = list(
    ids = uniprot_ids,
    from = "UniProtKB_AC-ID",
    to = "STRING"
  )
  
  # Make the initial submission
  r <- POST(url = "https://rest.uniprot.org/idmapping/run", body = files, encode = "multipart", accept_json())
  submission <- content(r, as = "parsed")
  print(paste("STRING mapping submission:", submission))
  
  if (isJobReady(submission[["jobId"]])) {
    Sys.sleep(10)
    url <- paste0("https://rest.uniprot.org/idmapping/results/", submission[["jobId"]])
    # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available
    url <- paste0(url, "?format=tsv")
    r <- GET(url = url)
    raw_content <- content(r, as = "raw")
    decompressed <- rawToChar(memDecompress(raw_content, type = "gzip"))
    stringTable <- read_tsv(I(decompressed))
    print(stringTable)
  }
  
  # Join STRING IDs back to original dataset
  nodes <- full_join(nodes, stringTable, by = c("uniprot" = "From"))
  nodes <- rename(nodes, string = To)
  head(nodes)
  
  # join STRING IDs back to edge list
  edges <- inner_join(edges, nodes, by = c("id1" = "id")) %>% inner_join(nodes, by = c("id2" = "id"))
  colnames(edges) <- c("id1", "id2", "uniprot1", "string1", "uniprot2", "string2")
  # extract only STRING IDs
  head(edges)
  string_edges <- edges %>% select(string1, string2)
  
  ### Get Ensembl IDs from Uniprot IDs
  ensembl_files = list(
    ids = uniprot_ids,
    from = "UniProtKB_AC-ID",
    to = "Ensembl"
  )
  
  # Make the initial submission
  r <- POST(url = "https://rest.uniprot.org/idmapping/run", body = ensembl_files, encode = "multipart", accept_json())
  submission <- content(r, as = "parsed")
  print(paste("ENSEMBL mapping submission:", submission))
  
  # When job is ready print the results
  if (isJobReady(submission[["jobId"]])) {
    Sys.sleep(10)
    url <- paste0("https://rest.uniprot.org/idmapping/results/", submission[["jobId"]])
    # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available
    url <- paste0(url, "?format=tsv")
    r <- GET(url = url)
    raw_content <- content(r, as = "raw")
    decompressed <- rawToChar(memDecompress(raw_content, type = "gzip"))
    ensemblTable <- read_tsv(I(decompressed))
    print(ensemblTable)
  }
  
  # Join STRING IDs back to original dataset
  nodes <- full_join(nodes, ensemblTable, by = c("uniprot" = "From"))
  nodes <- rename(nodes, ensembl = To)
  nodes$ensembl <- sapply(str_split(nodes$ensembl, "\\."), `[`, 1)
  
  
  write_tsv(string_edges, paste0("STRING-EDGES-", pathway_path))
  write_tsv(edges_rm, paste0("DEL-EDGES-", pathway_path))
  write_tsv(nodes, paste0("NODES-", pathway_path))
  write_tsv(nodes_rm, paste0("DEL-NODES-", pathway_path))
  
  # Function for taking a processed Node dataset from process_panther_pathway and checking for TFs
  process_TFs <- function(nodes, pathway_path) {
    humanTFs <- read_tsv('../human-interactome/Homo_sapiens_TF.txt')
    matches <- inner_join(nodes, humanTFs, by = c("ensembl" = "Ensembl"))
    sources <- matches %>% select(id, uniprot, string, ensembl)
    
    write_tsv(sources, paste0("SOURCES-", pathway_path))
    
    return(sources)
  }
  
  # Function for taking a processed Node dataset from process_panther_pathway and checking for receptors
  process_receptors <- function(nodes, pathway_path) {
    receptors <- read_tsv('../human-interactome/Homo_sapiens_surfaceome.txt')
    receptors <- receptors %>% 
      select(`UniProt accession`, `Ensembl gene`, `Membranome Almen main-class`) %>%
      filter(`Membranome Almen main-class` == "Receptors")
    matches <- inner_join(nodes, receptors, by = c("uniprot" = "UniProt accession"))
    targets <- matches %>% select(id, uniprot, string, ensembl)
    
    write_tsv(targets, paste0("TARGETS-", pathway_path))
    
    return(targets)
  }
  
  process_prizes <- function(sources, targets, pathway_path) {
    prizes100 <- rbind(sources, targets) %>% 
      mutate(prizes = 100, active = "true") %>% 
      select(id, prizes, active)
    
    write_tsv(prizes100, paste0("PRIZES-100-", pathway_path))
  }
  
  targets <- process_receptors(nodes, pathway_path)
  sources <- process_TFs(nodes, pathway_path)
  prizes <- process_prizes(sources, targets, pathway_path)
  
  # Return all four datasets as a list
  return(list(edges = string_edges, edges_rm = edges_rm, nodes = nodes, nodes_rm = nodes_rm, sources = sources, targets = targets, prizes = prizes))
}

# Example usage:
result <- process_panther_pathway(pathway_path)

# Access individual datasets
edges <- result$edges
edges_rm <- result$edges_rm
nodes <- result$nodes
nodes_rm <- result$nodes_rm
sources <- result$sources
targets <- result$targets
prizes <- result$prizes
