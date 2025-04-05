library(tidyverse)
# setwd('synthetic-data/')

# 13 million edges
human <- read_delim('human-interactome/9606.protein.links.full.v12.0.txt', delim = ' ')

head(human)

# 5 million edges
experiments_1 <- human %>% filter(experiments >= 1 | experiments_transferred >= 1)

# 3 million edges
experiments_100 <- human %>% filter(experiments >= 100 | experiments_transferred >= 100)

# 1 million edges
experiments_200 <- human %>% filter(experiments >= 200 | experiments_transferred >= 200)

# 500k
experiments_300 <- human %>% filter(experiments >= 300 | experiments_transferred >= 300)

# 300k
experiments_400 <- human %>% filter(experiments >= 400 | experiments_transferred >= 400)

# 200k
experiments_500 <- human %>% filter(experiments >= 500 | experiments_transferred >= 500)

# 200k
experiments_500 <- human %>% filter(experiments >= 500 | experiments_transferred >= 500)

experiments_600 <- human %>% filter(experiments >= 600 | experiments_transferred >= 600)

experiments_700 <- human %>% filter(experiments >= 700 | experiments_transferred >= 700)

experiments_800 <- human %>% filter(experiments >= 800 | experiments_transferred >= 800)

experiments_900 <- human %>% filter(experiments >= 900 | experiments_transferred >= 900)

hist(human$experiments)
hist(human$experiments_transferred)

write_tsv(experiments_1, 'interactomes/threshold-interactomes/human_interactome-1.txt')
write_tsv(experiments_100, 'interactomes/threshold-interactomes/human_interactome-100.txt')
write_tsv(experiments_200, 'interactomes/threshold-interactomes/human_interactome-200.txt')
write_tsv(experiments_300, 'interactomes/threshold-interactomes/human_interactome-300.txt')
write_tsv(experiments_400, 'interactomes/threshold-interactomes/human_interactome-400.txt')
write_tsv(experiments_500, 'interactomes/threshold-interactomes/human_interactome-500.txt')
write_tsv(experiments_600, 'interactomes/threshold-interactomes/human_interactome-600.txt')
write_tsv(experiments_700, 'interactomes/threshold-interactomes/human_interactome-700.txt')
write_tsv(experiments_800, 'interactomes/threshold-interactomes/human_interactome-800.txt')
write_tsv(experiments_900, 'interactomes/threshold-interactomes/human_interactome-900.txt')
