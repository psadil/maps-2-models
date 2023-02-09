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
source(here::here("R","hcp.R"))

Sys.setenv(TAR_PROJECT = "hcp")

targets::tar_load(prop0all, store="_hcp")
targets::tar_load(gray, store="_hcp")
targets::tar_load(at, store="_hcp")

tmp <- prop0all |>
  distinct(n_sub, Task) |>
  collect() |>
  expand(n_sub, Task)

pop_d <- targets::tar_read(pop_d, store="_hcp") |>
  semi_join(gray, by = c("x", "y", "z"), copy=TRUE) |>
  filter(!fct_match(d, "0"), sigma>0) |>
  mutate(d = cope / sigma) |>
  left_join(at, by=c("x", "y", "z")) |>
  select(-cope, -sigma, -n_sub, -ContrastName, -CopeNumber) |>
  filter(!is.na(`Network Name`)) |>
  right_join(tmp, by="Task") 

dd <- pop_d |>
  group_by(Task, `Network Name`, label) |>
  summarise(
    d = mean(d), 
    N = n(),
    .groups = "drop")

d2 <- prop0all  |>
  count(x, y, z, n_sub, Task) |>
  collect() |>
  right_join(pop_d, by = c("x", "y", "z", "n_sub", "Task")) |>
  filter(!is.na(`Network Name`)) |>
  mutate(sig = if_else(is.na(n), 0L, n) / 100) |>
  group_by(n_sub, Task, `Network Name`, label) |>
  summarise( # across voxels in each group, what proportion is over 80%
    s = sum(sig > 0.8),
    avg = mean(sig > 0.8), #s / n(),
    med = qbeta(.5, .5 + s - 1, .5 + n() - s - 1),
    lower = qbeta(.025, .5 + s - 1, .5 + n() - s - 1),
    upper = qbeta(.975, .5 + s - 1, .5 + n() - s - 1),
    .groups = "drop"
  ) |>
  mutate(
    lower=if_else(is.na(lower), 0, lower),
    upper=if_else(is.na(upper), 0, upper)) |>
  left_join(dd)

# 
# p1 <- d2 |>
#   mutate(
#     upper=if_else(is.na(lower), 0, lower),
#     lower=if_else(is.na(upper), 0, upper)) |>
#   ggplot(aes(x=n_sub, y=avg, color=`Network Name`)) +
#   facet_wrap(~Task) +
#   geom_line() +
#   geom_point() +
#   scale_x_continuous(
#     name = "N Sub",
#     breaks = seq(0, 100, by=20),
#     labels = seq(0, 100, by=20)
#   ) +
#   scale_y_continuous(
#     "proportion voxels > 80%",
#     breaks = c(0, .5, 1),
#     labels = c(0, .5, 1)) +
#   theme_gray(base_size = 9)
# 
# p2 <- pop_d |>
#   filter(!is.na(`Network Name`)) |>
#   count(Task, d, `Network Name`) |>
#   ggplot(aes(x=`Network Name`, y=n)) +
#   facet_wrap(~Task) +
#   geom_col(aes(color=`Network Name`, fill=factor(d)), position=position_dodge())
# 
# 
# p1 / p2 + plot_layout(guides = "collect")


d2 |>
  ggplot(aes(x=`Network Name`, y=avg, color=abs(d))) +
  facet_grid(Task~n_sub) +
  geom_point(aes(size=N)) +
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  scale_color_viridis_c(option="turbo") +
  scale_y_continuous(
    "proportion voxels > 80%",
    breaks = c(0, .5, 1),
    labels = c(0, .5, 1)) +
  theme_gray(base_size = 9) +
  theme(legend.position = "bottom")


p <- d2 |>
  ggplot(aes(x=abs(d), y=avg, color=`Network Name`)) +
  facet_grid(Task~n_sub) +
  geom_point(alpha=0.2) +
  # geom_errorbar(aes(ymin=lower, ymax=upper)) +
  scale_y_continuous(
    "proportion voxels > 80%",
    breaks = c(0, .5, 1),
    labels = c(0, .5, 1)) +
  theme(legend.position = "bottom")

ggsave(here::here("analyses/prop-active-by-network.png"), p, device = ragg::agg_png)


doit <- function(delta, n_voxels, N=60){
  power <- power.t.test(
    n = N, 
    delta=delta, 
    type="one.sample")
  qbinom(0.8, size=n_voxels, prob=power$power) / n_voxels
}

d2 |> 
  mutate(x = doit(d, n_voxels=N, N=n_sub)) |>
  ggplot(aes(x=d)) +
  geom_line(aes(y=x)) +
  geom_point(aes(y=med), alpha=0.2) +
  facet_wrap(~n_sub)
  
  
p <- d2 |>
  filter(n_sub==60) |>
  ggplot(aes(x=abs(d), y=avg)) +
  facet_grid(Task~`Network Name`) +
  geom_point(alpha=0.2) +
  xlab("|Cohen's d|") +
  scale_y_continuous(
    "proportion voxels > 80%",
    breaks = c(0, .5, 1),
    labels = c(0, .5, 1),
    limits = c(0, 1)) +
  theme(legend.position = "bottom")

ggsave(here::here("analyses/prop-active-by-roi.png"), p, device = ragg::agg_png)
