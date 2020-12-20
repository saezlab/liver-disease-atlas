library(tidyverse)
library(here)

# create annotations for the chronic and human contrasts ------------------
contrast_anno = tribble(
  ~contrast, ~disease, ~author, ~label, ~year, ~source, ~phenotype,
  "advanced_vs_mild", "NAFLD", "Diehl", "NAFLD (Diehl et al.)", NA_integer_, "diehl", "nafld",
  "hcv_adv_vs_early", "HCV", "Ramnath", "HCV (Ramnath et al.)",NA_integer_, "ramnath", "fibrosis",
  "nafld_adv_vs_early", "NAFLD", "Ramnath", "NAFLD (Ramnath et al.)",NA_integer_, "ramnath", "fibrosis",
  "stage_1_vs_0", "NAFLD Stage 1", "Hoang", "NAFLD Stage 1 (Hoang et al.)",NA_integer_, "hoang", "nafld",
  "stage_2_vs_0", "NAFLD Stage 2", "Hoang", "NAFLD Stage 2 (Hoang et al.)",NA_integer_, "hoang", "nafld",
  "stage_3_vs_0", "NAFLD Stage 3", "Hoang", "NAFLD Stage 3 (Hoang et al.)",NA_integer_, "hoang", "nafld",
  "stage_4_vs_0", "NAFLD Stage 4", "Hoang", "NAFLD Stage 4 (Hoang et al.)",NA_integer_, "hoang", "nafld",
  "stage_5_vs_0", "NAFLD Stage 5", "Hoang", "NAFLD Stage 5 (Hoang et al.)",NA_integer_, "hoang", "nafld",
  "stage_6_vs_0", "NAFLD Stage 6", "Hoang", "NAFLD Stage 6 (Hoang et al.)",NA_integer_, "hoang", "nafld",
  "steatosis_vs_ctrl", "NAFLD", "Hampe", "NAFLD (Hampe et al. '13)", 2013, "hampe13", "nash",
  "nash_vs_ctrl", "NASH", "Hampe", "NASH (Hampe et al. '13)", 2013, "hampe13", "nash",
  "nash_vs_ctrl", "NASH", "Hampe", "NASH (Hampe et al. '14)", 2014, "hampe14", "omni",
  "nafld_vs_ctrl", "NAFLD", "Hampe", "NAFLD (Hampe et al. '14)", 2014, "hampe14", "omni",
  "pbc_vs_ctrl", "PBC", "Hampe", "PBC (Hampe et al. '14)", 2014, "hampe14", "omni",
  "psc_vs_ctrl", "PSC", "Hampe", "PSC (Hampe et al. '14)", 2014, "hampe14", "omni"
) %>%
  mutate(label = fct_inorder(label)) %>%
  mutate(author2 = str_extract(label, "\\(.*\\)"),
         author2 = str_remove_all(author2, "\\(|\\)"))


saveRDS(contrast_anno,
        here("data/meta-mouse-vs-human/contrast_annotation.rds"))
