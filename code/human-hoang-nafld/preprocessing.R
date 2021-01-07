# Copyright (c) [2021] [Christian H. Holland]
# christian.holland@bioquant.uni-heidelberg.de

library(tidyverse)
library(here)
library(GEOquery)

source(here("code/utils-utils.R"))

# Download meta data ------------------------------------------------------
df <- getGEO("GSE130970")

meta_data <- df$GSE130970_series_matrix.txt.gz %>%
  pData() %>%
  rownames_to_column("sample") %>%
  as_tibble() %>%
  select(
    sample = title, fibrosis = `fibrosis stage:ch1`,
    lobular_inflammation = `lobular inflammation grade:ch1`,
    nafld = `nafld activity score:ch1`, gender = `Sex:ch1`,
    steatosis = `steatosis grade:ch1`, age = `age at biopsy:ch1`,
    ballooning = `cytological ballooning grade:ch1`
  ) %>%
  mutate( # sample = fct_inorder(sample),
    fibrosis = fct_inseq(fibrosis),
    lobular_inflammation = fct_inseq(lobular_inflammation),
    nafld = fct_inseq(nafld),
    gender = as_factor(str_to_lower(gender)),
    steatosis = fct_inseq(steatosis),
    age = as.numeric(age),
    ballooning = fct_inseq(ballooning)
  ) %>%
  mutate(nafld = factor(str_c("stage", nafld, sep = "_"),
    levels = str_c("stage", c(0:6), sep = "_")
  ))



# Load count matrix -------------------------------------------------------
# read count data with entrez gene ids
count_mat_entrez <- read_csv(
  here("data/human-hoang-nafld/GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv")
)

# translate gene entrez ids to hgnc symbols
count_mat_hgnc <- count_mat_entrez %>%
  mutate(entrez_id = as.character(entrez_id)) %>%
  rename(gene = entrez_id) %>%
  translate_gene_ids(from = "entrez_hgnc", to = "symbol_hgnc")

# summarize counts for duplicated genes (takes quite long)
count_mat <- count_mat_hgnc %>%
  drop_na(gene) %>%
  group_by(gene) %>%
  summarise(across(where(is.numeric), median)) %>%
  data.frame(row.names = 1, check.names = F)

# save meta data and count matrix -----------------------------------------
stopifnot(meta_data$sample == colnames(count_mat))

saveRDS(meta_data, here("data/human-hoang-nafld/meta_data.rds"))
saveRDS(count_mat, here("data/human-hoang-nafld/count_matrix.rds"))
