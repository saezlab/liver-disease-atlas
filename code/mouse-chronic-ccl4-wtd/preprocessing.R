# Copyright (c) [2021] [Christian H. Holland]
# christian.holland@bioquant.uni-heidelberg.de

library(tidyverse)
library(readxl)
library(here)

# load count matrix from txt file -----------------------------------------
count_matrix <- read_delim(
  here("data/mouse-chronic-ccl4-wtd/GSE99010-GPL17021_series_matrix.txt"),
  delim = "\t", comment = "!"
) %>%
  # remove samples related to cancer
  select(-c("GSM2630033", "GSM2630034")) %>%
  data.frame(row.names = 1, check.names = F)

# manually build meta data with information from GEO ----------------------
meta_data <- tibble(sample = colnames(count_matrix)) %>%
  mutate(
    diet = as_factor(c("nd", "nd", "wd", "wd", "nd", "nd", "wd", "wd")),
    treatment = as_factor(c(
      "oil", "ccl4", "oil", "ccl4", "oil", "ccl4", "oil", "ccl4"
    )),
    time = ordered(c(12, 12, 12, 12, 24, 24, 24, 24))
  ) %>%
  unite(group, diet, treatment, time, remove = FALSE) %>%
  mutate(group = as_factor(group))

# save meta data and count matrix -----------------------------------------
stopifnot(meta_data$sample == colnames(count_matrix))

saveRDS(meta_data, here("data/mouse-chronic-ccl4-wtd/meta_data.rds"))
saveRDS(count_matrix, here("data/mouse-chronic-ccl4-wtd/count_matrix.rds"))
