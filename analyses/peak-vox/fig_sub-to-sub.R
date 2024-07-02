library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

Sys.setenv(TAR_PROJECT = "hcp")

# what proportion of studies have a peak that is within x mm of the gold standard?


gold_peaks <- targets::tar_read(gold_peaks) |>
  select(Task, m) |>
  mutate(m = map(m, ~ mutate(.x, peak_i = row_number(desc(Value))))) |>
  unnest(m) |>
  filter(peak_i < 11) |>
  distinct(Task, x, y, z, peak_i)

# TODO: get closest peaks in all individual images ...
