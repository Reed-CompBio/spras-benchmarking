library(dplyr)
library(tidyr)

# match path with where the pathway is at
pathway_path <- "/Users/abarelvi/spras-benchmarking/string/Beta3_adrenergic/Beta3_adrenergic_rec.txt"

df_wnt <- read.csv(pathway_path, header = TRUE, sep ="\t")
df_wnt <- df_wnt %>%
  select(PARTICIPANT_A, PARTICIPANT_B) %>%
  rename(id1 = PARTICIPANT_A) %>%
  rename(id2 = PARTICIPANT_B)

genes <- data.frame(unique(c(df_wnt$id1, df_wnt$id2)))

write.csv(genes, "/Users/abarelvi/spras-benchmarking/string/genes.csv", row.names = FALSE)

# use this output file to map from Gene to UniProtKB/Swiss-Prot on human through https://www.uniprot.org/id-mapping
# name is file gene_to_uniprot.tsv

df_map_gene_uni <- read.csv("/Users/abarelvi/spras-benchmarking/string/Beta3_adrenergic/gene_to_uniprot.tsv", header = TRUE, sep = "\t")
df_map_gene_uni <- df_map_gene_uni %>% 
  select(From, Entry) %>%
  rename( gene = From) %>%
  rename( uniprot = Entry)

uniprot <- data.frame(unique(c(df_map_gene_uni$uniprot)))

write.csv(uniprot,  "/Users/abarelvi/spras-benchmarking/string/uniprot.csv", row.names = FALSE)

# use this output file to map from UniProtKB/Swiss-Prot to STRING human through https://www.uniprot.org/id-mapping
# name is file uniprot_to_string.tsv

df_map_uni_string <- read.csv("/Users/abarelvi/spras-benchmarking/string/Beta3_adrenergic/uniprot_to_string.tsv", header = TRUE, sep = "\t")
df_map_uni_string <- df_map_uni_string %>%
  rename(uniprot = From) %>%
  rename(string = To)

set.seed(123)  # Set seed for reproducibility
df_map_gene_uni_one <- df_map_gene_uni %>%
  group_by(gene) %>%
  slice_sample(n = 1) %>%  # Randomly select one row per gene
  ungroup()

many_to_many_gene <- df_map_gene_uni_one %>%
  group_by(gene) %>%
  summarise(count_uniprot = n_distinct(uniprot)) %>%
  filter(count_uniprot > 1)


set.seed(123)  # Set seed for reproducibility
df_map_uni_string_one <- df_map_uni_string %>%
  group_by(uniprot) %>%
  slice_sample(n = 1) %>%  # Randomly select one row per gene
  ungroup()

many_to_many_string <- df_map_uni_string_one %>%
  group_by(uniprot) %>%
  summarise(count_string = n_distinct(string)) %>%
  filter(count_string > 1)


df_wnt_mapped <- df_wnt %>%
  left_join(df_map_gene_uni_one, by = c("id1" = "gene")) %>%
  rename(uniprot_id1 = uniprot) %>%
  left_join(df_map_gene_uni_one, by = c("id2" = "gene")) %>%
  rename(uniprot_id2 = uniprot)

# Step 2: Left join the result with df_map_uni_string to get String IDs for uniprot_id1 and uniprot_id2
df_final <- df_wnt_mapped %>%
  left_join(df_map_uni_string_one, by = c("uniprot_id1" = "uniprot")) %>%
  rename(string_id1 = string) %>%
  left_join(df_map_uni_string_one, by = c("uniprot_id2" = "uniprot")) %>%
  rename(string_id2 = string)

df_final <- df_final %>%
  drop_na()

df_final_clean <- df_final[!duplicated(df_final), ]


write.csv(df_final_clean,  "/Users/abarelvi/spras-benchmarking/string/wnt_mapped.csv", row.names = FALSE )
