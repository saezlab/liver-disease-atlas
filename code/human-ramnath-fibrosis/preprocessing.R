# Copyright (c) [2021] [Christian H. Holland]
# christian.holland@bioquant.uni-heidelberg.de

library(tidyverse)
library(here)

# Download meta data ------------------------------------------------------
meta_raw <- read_delim("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6863/E-MTAB-6863.sdrf.txt", delim = "\t")

meta_data <- meta_raw %>%
  select(
    sample = "Source Name", age = "Characteristics[age]",
    gender = "Characteristics[sex]", disease = "Characteristics[disease]",
    stage = "Characteristics[disease staging]",
    infect = "Characteristics[infect]",
    bmi = "Characteristics[body mass index]",
    steatosis = "Characteristics[clinical information]",
    fibrosis = "Characteristics[fibrosis stage]"
  ) %>%
  mutate(disease = case_when(
    disease == "liver disease" ~ "hcv",
    disease == "non-alcoholic fatty liver disease" ~ "nafld",
    disease == "alcoholic liver disease" ~ "ald",
    TRUE ~ disease
  )) %>%
  mutate(
    stage = as_factor(stage),
    gender = as_factor(gender),
    disease = as_factor(disease),
    infect = as_factor(infect),
    # raises a warning
    bmi = as.numeric(bmi),
    steatosis = as_factor(steatosis),
    fibrosis = as_factor(fibrosis)
  ) %>%
  mutate(
    group = str_c(disease, stage, sep = "_"),
    group = as_factor(group)
  ) %>%
  mutate(sample_num = as.character(parse_number(sample))) %>%
  mutate(sample_num = case_when(
    nchar(sample_num) == 1 ~ str_pad(sample_num, 2,
      pad = 0
    ),
    TRUE ~ sample_num
  )) %>%
  mutate(
    sample = str_c("sample", sample_num),
    sample_num = NULL
  ) %>%
  arrange(sample) %>%
  # remove patients with alcohol liver disease
  filter(disease != "ald") %>%
  mutate(
    disease = fct_drop(disease),
    group = fct_drop(group)
  )

# Consolidate count matrix ------------------------------------------------
count_matrix_all <- list.files(here("data/human-ramnath-fibrosis"),
  full.names = T, pattern = ".txt"
) %>%
  map(function(path) {
    sample_id <- str_extract(path, "sample[0-9]*")
    read_delim(path, delim = "\t", col_names = c("gene", sample_id))
  }) %>%
  reduce(left_join) %>%
  data.frame(row.names = 1, check.names = F)

count_matrix <- count_matrix_all[, meta_data$sample]

# save meta data and count matrix -----------------------------------------
stopifnot(meta_data$sample == colnames(count_matrix))

saveRDS(meta_data, here("data/human-ramnath-fibrosis/meta_data.rds"))
saveRDS(count_matrix, here("data/human-ramnath-fibrosis/count_matrix.rds"))
