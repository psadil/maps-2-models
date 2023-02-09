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

space <- targets::tar_read(space, store=here::here("_hcp"))

space |> 
  na.omit() |>
  filter(corrp_thresh%in%c(0.95)) |> 
  mutate(
    Value = Value / sqrt(n_pop),
    Value = cut(Value, breaks=c(0, 0.2, 0.5, 0.8, Inf), right=FALSE)) |>
  whoppeR::WISEsummary(
    dependentvars = "d",
    betweenvars = c("n_sub", "Value", "corrp_thresh"),
    na.rm = TRUE
  ) |>
  mutate(CI=map2_chr(d_CI_lower, d_CI_upper, ~str_c(round(.x,digits=2), round(.y,digits=2), sep = ", "))) |>
  select(N=n_sub, FWE=corrp_thresh, `Effect Size`=Value, `Average Distance`=d_mean, SEM=d_sem, Voxels=d_n, CI) |>
  print(n=100) |>
  knitr::kable(format = "latex", digits = 2)


gray <- to_tbl(MNITemplate::getMNISegPath(res="2mm")) |>
  dplyr::filter(value==2)

targets::tar_load(pop_d, store="_hcp")

pop_d |>
  filter(fct_match(d, c("small", "medium", "large"))) |>
  group_by(ContrastName) |>
  mutate(N = n()) |>
  group_by(ContrastName, d) |>
  summarise(
    prop = n() / nrow(gray) * 100,
    .groups = "drop"
  ) |>
  pivot_wider(names_from = ContrastName, values_from = prop) |>
  knitr::kable(format = "latex", digits = 2)


p <- pop_d |> 
  filter(str_detect(ContrastName, "PUNISH-REWARD", negate = TRUE)) |>
  # na.omit() |>
  # filter(!fct_match(d, "0")) |>
  ggplot() + 
  facet_wrap(~Task) +
  # geom_point(aes(x=cope, y=sigma), alpha=0.1) +
  scattermore::geom_scattermore(aes(x=cope, y=sigma), alpha=0.1, pointsize=1) +
  geom_function(fun = function(x) sqrt(20)*x/qt(0.001, 19), xlim=c(-200,0)) + 
  geom_function(fun = function(x) sqrt(100)*x/qt(0.001, 99), xlim=c(-200,0), color = "gray50") + 
  geom_function(fun = function(x) sqrt(100)*x/qt(0.999, 99), xlim=c(0,200), color = "gray50") + 
  geom_function(fun = function(x) sqrt(20)*x/qt(0.999, 19), xlim=c(0,200)) + 
  coord_cartesian(xlim=c(-200, 200), ylim = c(0, 150)) +
  xlab(expression(beta~mean))  +
  ylab(expression(beta~SD)) +
  theme_gray(base_size = 10)
ggsave("analyses/beta-by-sd.png", device = ragg::agg_png, width = 6, height = 4)


targets::tar_load(prop0all, store="_hcp")

tmp <- prop0all |>
  distinct(n_sub, Task) |>
  collect() |>
  expand(n_sub, Task)

pop_d <- targets::tar_read(pop_d, store="_hcp") |>
  semi_join(gray, by = c("x", "y", "z")) |>
  filter(!fct_match(d, "0"), sigma>0) |>
  select(-cope, -sigma, -n_sub, -ContrastName, -CopeNumber) |>
  right_join(tmp, by="Task")

d2 <- targets::tar_read(prop0all, store="_hcp")  |>
  count(x, y, z, n_sub, Task) |>
  collect() |>
  right_join(pop_d, by = c("x", "y", "z", "n_sub", "Task")) |>
  mutate(sig = if_else(is.na(n), 0L, n) / 100) |>
  group_by(n_sub, d, Task) |>
  summarise( # across voxels in each group, what proportion is over 80%
    s = sum(sig > 0.8),
    avg = mean(sig > 0.8), #s / n(),
    med = qbeta(.5, .5 + s - 1, .5 + n() - s - 1),
    lower = qbeta(.025, .5 + s - 1, .5 + n() - s - 1),
    upper = qbeta(.975, .5 + s - 1, .5 + n() - s - 1),
    .groups = "drop"
  ) 

p2 <- d2 |>
  ggplot(aes(x=n_sub, y=avg, color=d)) +
  facet_wrap(~Task) +
  geom_point() +
  geom_line() +
  guides(color = guide_legend(title="Effect Size")) +
  scale_x_continuous(
    name = "N Sub",
    breaks = seq(0, 100, by=20),
    labels = seq(0, 100, by=20)
  ) +
  ylab("proportion voxels > 80%") +
  theme_gray(base_size = 9)

ggsave("analyses/prop-active.png", p2, width=6, height = 2)


d3 <- targets::tar_read(prop0all, store="_hcp")  |>
  count(x, y, z, n_sub, Task) |>
  collect() |>
  right_join(pop_d, by = c("x", "y", "z", "n_sub", "Task")) |>
  mutate(sig = 1 - if_else(is.na(n), 0L, n) / 100) |>
  group_by(n_sub, d, Task) |>
  summarise( 
    avg = mean(sig),
    .groups = "drop"
  ) 

p3 <- d3 |>
  ggplot(aes(x=n_sub, y=avg, color=d)) +
  facet_wrap(~Task) +
  geom_point() +
  geom_line() +
  scale_x_continuous(
    name = "N Sub",
    breaks = seq(0, 100, by=20),
    labels = seq(0, 100, by=20)
  )  +
  guides(color = guide_legend(title="Effect Size")) +
  ylab("false negative rate") +
  theme_gray(base_size = 9)

ggsave("analyses/false-neg_hcp.png", p3, width=6, height = 2)
