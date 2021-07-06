# Copyright (c) [2021] [Christian H. Holland]
# christian.holland@bioquant.uni-heidelberg.de

library(tidyverse)
library(readxl)
library(here)

# load count matrix from txt file -----------------------------------------
count_matrix <- list.files(
  here("data/mouse-chronic-ccl4-etoh"),
  pattern = "txt.gz", full.names = TRUE
) %>%
  # remove samples related to cancer
  discard(str_detect(
    .,
    "GSM3375274|GSM3375275|GSM3375276|GSM3375277|GSM3375278"
  )) %>%
  map_dfc(function(path) {
    sample <- str_extract(path, "GSM\\d+")
    read_delim(path, delim = "\t", col_names = c("gene", sample)) %>%
      data.frame(row.names = 1, check.names = F)
  })

# manually build meta data with information from GEO ----------------------
meta_data <- tibble(sample = colnames(count_matrix)) %>%
  mutate(
    time = ordered(c(0, 0, 20, 20, 20, 20, 20, 20, 20, 20)),
    treatment = as_factor(c(
      "control", "control", "etoh", "etoh", "etoh", "ccl4", "ccl4", "ccl4_etoh",
      "ccl4_etoh", "ccl4_etoh"
    ))
  ) %>%
  unite(group, treatment, time, remove = FALSE) %>%
  mutate(group = as_factor(group)) %>%
  select(sample, time, treatment, group)

# save meta data and count matrix -----------------------------------------
stopifnot(meta_data$sample == colnames(count_matrix))

saveRDS(meta_data, here("data/mouse-chronic-ccl4-etoh/meta_data.rds"))
saveRDS(count_matrix, here("data/mouse-chronic-ccl4-etoh/count_matrix.rds"))
