library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

Sys.setenv(TAR_PROJECT = "hcp")

targets::tar_load(pop_d)

pop_d |> 
  ggplot() + 
  facet_wrap(~Task) +
  scattermore::geom_scattermore(aes(x=cope, y=sigma), alpha=0.1, pointsize=1) +
  geom_function(fun = function(x) sqrt(20)*x/qt(0.001, 19), xlim=c(-200,0)) + 
  geom_function(fun = function(x) sqrt(100)*x/qt(0.001, 99), xlim=c(-200,0), color = "gray50") + 
  geom_function(fun = function(x) sqrt(100)*x/qt(0.999, 99), xlim=c(0,200), color = "gray50") + 
  geom_function(fun = function(x) sqrt(20)*x/qt(0.999, 19), xlim=c(0,200)) + 
  coord_cartesian(xlim=c(-200, 200), ylim = c(0, 150)) +
  xlab(expression(beta~mean))  +
  ylab(expression(beta~SD)) +
  theme_gray(base_size = 14)

ggsave("analyses/figures/topo_gold.png", device = ragg::agg_png, width = 8, height = 6)


pop_d |> 
  mutate(d = abs(cope / sigma)) |>
  na.omit() |>
  ggplot(aes(x=d, y=Task)) + 
  ggdist::stat_dots(quantiles = 1000) +
  xlab("Abs. Effect Size") 

ggsave(
  "analyses/figures/topo_abs-vox-eff.png", 
  device = ragg::agg_png, 
  width = 4, 
  height = 4)

