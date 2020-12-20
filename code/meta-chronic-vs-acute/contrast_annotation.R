library(tidyverse)
library(here)

# create annotations for the chronic and acute contrasts ------------------
contrast_annotation = tribble(
  ~contrast, ~treatment, ~time_label, ~time_label2, ~source, ~class,
  "treat_vs_ctrl", "Tunicamycin", NA_character_, NA_character_, "Wek", "Acute",
  "inLiver_lps_vs_ctrl", "LPS", NA_character_, NA_character_, "Godoy", "Acute",
  "ccl_2h_vs_0h", "CCl4", "2 Hours", "Hour 2", "Godoy", "Acute",
  "ccl_8h_vs_0h", "CCl4", "8 Hours", "Hour 8", "Godoy", "Acute",
  "ccl_24h_vs_0h", "CCl4", "1 Day", "Day 1", "Godoy", "Acute",
  "ccl_48h_vs_0h", "CCl4", "2 Days", "Day 2", "Godoy", "Acute",
  "ccl_96h_vs_0h", "CCl4", "4 Days", "Day 4", "Godoy", "Acute",
  "ccl_144h_vs_0h", "CCl4", "6 Days", "Day 6", "Godoy", "Acute",
  "ccl_192h_vs_0h", "CCl4", "8 Days", "Day 8", "Godoy", "Acute",
  "ccl_384h_vs_0h", "CCl4", "16 Days", "Day 16", "Godoy", "Acute",
  "apap_1h_vs_0h", "APAP", "1 Hour", "Hour 1", "Ghallab", "Acute",
  "apap_6h_vs_0h", "APAP", "6 Hours", "Hour 6", "Ghallab", "Acute",
  "apap_12h_vs_0h", "APAP", "12 Hours", "Hour 12", "Ghallab", "Acute",
  "apap_24h_vs_0h", "APAP", "1 Day", "Day 1", "Ghallab", "Acute",
  "apap_48h_vs_0h", "APAP", "2 Days", "Day 2", "Ghallab", "Acute",
  "apap_96h_vs_0h", "APAP", "4 Days", "Day 4", "Ghallab", "Acute",
  "apap_144h_vs_0h", "APAP", "6 Days", "Day 6", "Ghallab", "Acute",
  "apap_192h_vs_0h", "APAP", "8 Days", "Day 8", "Ghallab", "Acute",
  "apap_384h_vs_0h", "APAP", "16 Days", "Day 16", "Ghallab", "Acute",
  "ph_0.043d", "Partial Hepatectomy", "1 Hour", "Hour 1", "Godoy", "Acute",
  "ph_0.25d", "Partial Hepatectomy", "6 Hours", "Hour 6", "Godoy", "Acute",
  "ph_0.5d", "Partial Hepatectomy", "12 Hours", "Hour 12", "Godoy", "Acute",
  "ph_1d", "Partial Hepatectomy", "1 Day", "Day 1", "Godoy", "Acute",
  "ph_2d", "Partial Hepatectomy", "2 Days", "Day 2", "Godoy", "Acute",
  "ph_3d", "Partial Hepatectomy","3 Days", "Day 3", "Godoy", "Acute",
  "ph_4d", "Partial Hepatectomy","4 Days", "Day 4", "Godoy", "Acute",
  "ph_7d", "Partial Hepatectomy","1 Week", "Week 1", "Godoy", "Acute",
  "ph_14d", "Partial Hepatectomy","2 Weeks", "Week 2", "Godoy", "Acute",
  "ph_28d", "Partial Hepatectomy","1 Month", "Month 1", "Godoy", "Acute",
  "ph_84d", "Partial Hepatectomy","3 Months", "Month 3", "Godoy", "Acute",
  "bdl_vs_sham_1d", "Bile Duct Ligation", "1 Day", "Day 1", "Ghallab", "Acute",
  "bdl_vs_sham_3d","Bile Duct Ligation", "3 Days", "Day 3", "Ghallab", "Acute",
  "bdl_vs_sham_7d","Bile Duct Ligation", "7 Days", "Day 7", "Ghallab", "Acute",
  "bdl_vs_sham_21d","Bile Duct Ligation", "21 Days", "Day 21", "Ghallab", "Acute",
  "pure_ccl_2m_vs_0m", "CCl4", "2 Months", "Month 2", "Ghallab", "Chronic",
  "pure_ccl_6m_vs_0m", "CCl4", "6 Months", "Month 6", "Ghallab", "Chronic",
  "pure_ccl_12m_vs_0m", "CCl4", "12 Months", "Month 12", "Ghallab", "Chronic"
) %>%
  mutate(treatment_abbr = case_when(
    str_detect(treatment, "Bile|Partial") ~ abbreviate(treatment, 1),
    TRUE ~ treatment
  )) %>%
  mutate(time = str_match(contrast, "[0-9]*\\.*[0-9]+[a-z]{1}")[,1]) %>%
  # get time value and unit of timepoint
  mutate(value = parse_number(time),
         unit = str_extract(time, "[a-z]")) %>%
  # convert unit into days for ph experiment
  mutate(value = case_when(
    str_detect(contrast, "ph") ~ as.integer(value * 24),
    TRUE ~ as.integer(value)
  )) %>%
  mutate(unit = case_when(
    str_detect(contrast, "ph") ~ "h",
    TRUE ~ unit
  )) %>%
  # convert time into a period
  mutate(time = case_when(
    unit == "h" ~ lubridate::dhours(value),
    unit == "d" ~ lubridate::ddays(value),
    unit == "m" ~ lubridate::duration(value, "month"))) %>%
  mutate(time_label = fct_reorder(time_label, time),
         time_label2 = fct_reorder(time_label2, time),
         contrast = fct_inorder(contrast)) %>%
  mutate(label = case_when(
    is.na(time_label) ~ treatment,
    TRUE ~ str_c(treatment, " (", time_label, ")")
  )) %>%
  mutate(label2 = case_when(
    is.na(time_label) ~ treatment,
    TRUE ~ str_c(treatment_abbr, " (", time_label, ")"),
  )) %>%
  mutate(label = fct_inorder(label),
         label2 = fct_inorder(label2),
         treatment_abbr = fct_inorder(treatment_abbr))

saveRDS(contrast_annotation,
        here("data/meta-chronic-vs-acute/contrast_annotation.rds"))
