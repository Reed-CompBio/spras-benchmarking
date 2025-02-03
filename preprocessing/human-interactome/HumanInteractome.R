library(tidyverse)
setwd('~/Desktop/GitHub_Repos/spras-benchmarking/preprocessing/human-interactome/')

# 13 million edges
human <- read_delim('9606.protein.links.full.v12.0.txt', delim = ' ')

head(human)

# 5 million edges
experiments <- human %>% filter(experiments > 0 | experiments_transferred > 0)


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

hist(human$experiments)
hist(human$experiments_transferred)

write_tsv(experiments_500, 'human_interactome-500.txt')
