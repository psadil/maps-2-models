library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

# NOT RUN!
# have not yet configured analysis correctly
# need to run analysis with training on UKB (to get large enough sample size)
# .  but then testing on the same HCP population (for comparability)

Sys.setenv(TAR_PROJECT = "hcp_ptfce")

targets::tar_load(c(data_model_gold_gold_to_study, data_model_gold_gold_to_study_ukb))

ukb <- data_model_gold_gold_to_study_ukb |>
  filter(measure == "age") |>
  mutate(measure = case_match(measure, "age" ~ "Age_in_Yrs")) |>
  mutate(dataset = "UKB")

data_model_gold_gold_to_study |>
  filter(measure == "Age_in_Yrs", task == "EMOTION") |>
  mutate(dataset = "HCP") |>
  bind_rows(ukb) |>
  filter(type == "simulation") |>
  ggplot(aes(x = n_sub, y = avg, color = dataset)) +
  facet_wrap(~task) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  scale_x_log10()
