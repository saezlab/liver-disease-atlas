library(tidyverse)
library(readxl)
library(here)

# tidy meta data ----------------------------------------------------------
meta_data <- read_excel(
  here("data/mouse-chronic-ccl4/RNAseq_Kiel_samples_Ghallab.xlsx")
) %>%
  rename(id = "name,barcode,plate,row,col", treatment = dose) %>%
  mutate(
    sample = str_trunc(id, 6, ellipsis = ""),
    time = as.ordered(month),
    treatment = as_factor(treatment),
    group = paste0(treatment, ".", time),
    group = case_when(
      group == "ctrl.0" ~ "wt",
      TRUE ~ group
    )
  ) %>%
  select(sample, time, treatment, group) %>%
  mutate(group = factor(group, levels = c(
    "wt",
    "oil.2", "ccl4.2",
    "ccl4.6",
    "oil.12", "ccl4.12"
  )))

# load count matrix from txt file -----------------------------------------
count_matrix <- read_delim(
  here("data/mouse-chronic-ccl4/raw_read_counts.txt"),
  delim = "\t"
) %>%
  data.frame(row.names = 1, check.names = F)

# save meta data and count matrix -----------------------------------------
stopifnot(meta_data$sample == colnames(count_matrix))

saveRDS(meta_data, here("data/mouse-chronic-ccl4/meta_data.rds"))
saveRDS(count_matrix, here("data/mouse-chronic-ccl4/count_matrix.rds"))
