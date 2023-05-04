library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)
source(here::here("R","utils.R"))
source(here::here("R","tfce.R"))
source(here::here("R","updates.R"))
source(here::here("R","poster.R"))

Sys.setenv(TAR_PROJECT = "hcp")


targets::tar_read(pop_d) |>
  filter(!is.na(d)) |>
  count(Task, d, name = "N") |>
  group_by(Task) |>
  mutate(Proportion = N / sum(N)) |>
  ggplot(aes(x=Task, fill=d, y = Proportion)) +
  geom_col(position = "dodge") +
  guides(fill = guide_legend("Cohen's d")) +
  xlab(NULL) +
  theme_gray(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 

ggsave("analyses/figures/prop-effect-size.png", width=4.5, height = 3)
