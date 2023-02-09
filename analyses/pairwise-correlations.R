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

targets::tar_load(pop_cor)

# correlation with copes
p <- pop_cor |>
  mutate(`N Sub` = factor(n_sub)) |>
  ggplot(aes(y=rho, x=`N Sub`, color = method)) +
  facet_wrap(~Task) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha=0.25, position = position_jitterdodge()) +
  ylab("Correlation with Reference")

ggsave("analyses/cor_w_reference.png", p)

# this is correlation of effect size (cohen's d)
targets::tar_load(pop_cor2)
pop_cor2 |>
  mutate(`N Sub` = factor(n_sub)) |>
  ggplot(aes(y=rho, x=`N Sub`, color = method)) +
  facet_wrap(~Task) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha=0.25, position = position_jitterdodge()) +
  ylab("Correlation with Reference")

targets::tar_load(pairwise)

tmp <- pairwise |>
  filter(y > x) |>
  mutate(f = atanh(r)) |>
  group_by(Task, n_sub, ContrastName, method) |>
  summarise(
    f = mean(f),
    r = mean(r),
    N = n(),
    .groups = "drop"
  ) |>
  mutate(
    `N Sub` = factor(n_sub),
    SE = 1/sqrt(N-3)) 

p <- tmp |>
  mutate(
    rr = tanh(f),
    lower = tanh(f - 1.96*SE),
    upper = tanh(f + 1.96*SE)) |>
  ggplot(aes(x=rr, y=`N Sub`, color=Task)) +
  facet_wrap(~method) +
  geom_point() +
  geom_errorbarh(aes(xmin=lower, xmax=upper)) +
  xlab("Average Pairwise Correlation [95% CI]")

ggsave("analyses/pairwise-correlations.png", p, device = ragg::agg_png)


p <- pairwise |>
  filter(method=="spearman") |>
  mutate(
    y = as.numeric(y),
    x = as.numeric(x)) |>
  filter(y < x) |>
  mutate(
    n_sub = factor(n_sub),
    n_sub = forcats::fct_relabel(n_sub, .fun = ~glue::glue("N Sub: {.x}"))) |>
  ggplot(aes(x = y, y = x)) +
  facet_grid(n_sub~Task) +
  geom_raster(aes(fill=r)) +
  scale_fill_viridis_c(option="turbo") +
  ylab("Iteration") +
  guides(fill = guide_colourbar(title = "Correlation")) +
  scale_x_continuous(
    name = "Iteration",
    breaks = c(1, 100),
    labels = c(1, 100)
  ) +
  scale_y_reverse(
    name = "Iteration",
    breaks = c(1, 100),
    labels = c(1, 100)
  ) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  )

ggsave("analyses/correlation-matrices.png", p, device = ragg::agg_png)



make_sd_plot <- function(gray, task){
  tmp <- gray |>
    filter((z %% 5) == 0, between(z, 15, 70))  |>
    filter(.data$Task == .env$task) |>
    collect() |>
    mutate(d = cope / sigma * correct_d(n_sub)) |>
    group_by(n_sub, ContrastName, Task, x, y, z) |>
    summarise(
      s = sd(d),
      .groups = "drop"
    ) 
  
  tmp |>
    mutate(`N Sub` = n_sub) |>
    ggplot(aes(x=x, y=y)) +
    coord_fixed() +
    facet_grid(`N Sub` ~ z, labeller = label_both) +
    geom_raster(
      aes(fill=value),
      data = 
        to_tbl(MNITemplate::getMNIPath(what="Brain", res = "2mm")) |>
        filter(z %in% unique(tmp$z)), show.legend = FALSE) +
    scale_fill_distiller(
      type = "seq",
      direction = -1,
      palette = "Greys") +
    ggnewscale::new_scale_fill() +
    geom_raster(aes(fill=s), alpha = 0.5) +
    scale_fill_viridis_c(
      option="turbo",
      limits = c(0, 0.5)) +
    scale_alpha_continuous(guide="none",range = c(0.3,1)) +
    labs(fill="sd(Study g)") +
    theme_void() +
    ggtitle(task)
}

ps <- map(unique(tfce$Task), ~make_sd_plot(gray, .x))

for (x in ps){
  ggsave(
    glue::glue("analyses/{x$labels$title}.png"),
    x,
    device = ragg::agg_png,
    dpi = 600)
}


make_avg_plot <- function(gray, task){
  tmp <- gray |>
    filter((z %% 5) == 0, between(z, 15, 70)) |>
    filter(.data$Task == .env$task) |>
    collect() |>
    mutate(d = cope / sigma * correct_d(n_sub)) |>
    group_by(n_sub, ContrastName, Task, x, y, z) |>
    summarise(
      s = mean(d, na.rm=TRUE),
      .groups = "drop"
    ) 
  
  tmp |>
    mutate(`N Sub` = n_sub) |>
    ggplot(aes(x=x, y=y)) +
    coord_fixed() +
    facet_grid(`N Sub` ~ z, labeller = label_both) +
    geom_raster(
      aes(fill=value),
      data = 
        to_tbl(MNITemplate::getMNIPath(what="Brain", res = "2mm")) |>
        filter(z %in% unique(tmp$z)), show.legend = FALSE) +
    scale_fill_distiller(
      type = "seq",
      direction = -1,
      palette = "Greys") +
    ggnewscale::new_scale_fill() +
    geom_raster(aes(fill=s), alpha = 0.5) +
    scale_fill_viridis_c(
      option="turbo",
      limits = c(-2, 2)) +
    scale_alpha_continuous(guide="none",range = c(0.3,1)) +
    labs(fill="avg(Study g)") +
    theme_void() +
    ggtitle(task)
}

pa <- map(unique(tfce$Task), ~make_avg_plot(gray, .x))

for (x in pa){
  ggsave(
    glue::glue("analyses/{x$labels$title}_avg.png"),
    x,
    device = ragg::agg_png,
    dpi = 600)
}


# plot d by variance
crossing(
  N = unique(tfce$n_sub),
  d = seq(0, 2, length.out=100)) |>
  mutate(
    v = d_var(N, d),
    N = factor(N),
  ) |>
  ggplot(aes(x=d, y=v)) + 
  geom_line(aes(color = N))


# differences
tmp <- gray |>
  filter((z %% 5) == 0, between(z, 15, 70), n_sub %in% c(20, 100)) |>
  collect() |>
  mutate(d = cope / sigma * correct_d(n_sub)) |>
  group_by(n_sub, ContrastName, Task, x, y, z) |>
  summarise(
    s = mean(d, na.rm=TRUE),
    .groups = "drop"
  ) |> 
  pivot_wider(names_from = n_sub, values_from = s) |> 
  mutate(d = `100` - `20`) 

p <- tmp |>
  ggplot(aes(x=x, y=y)) +
  coord_fixed() +
  facet_grid(Task ~ z, labeller = label_both) +
  geom_raster(
    aes(fill=value),
    data = 
      to_tbl(MNITemplate::getMNIPath(what="Brain", res = "2mm")) |>
      filter(z %in% unique(tmp$z)), show.legend = FALSE) +
  scale_fill_distiller(
    type = "seq",
    direction = -1,
    palette = "Greys") +
  ggnewscale::new_scale_fill() +
  geom_raster(aes(fill=d), alpha = 0.5) +
  scale_fill_viridis_c(
    option="turbo") +
  scale_alpha_continuous(guide="none",range = c(0.3,1)) +
  labs(fill="g_100 - g_20") +
  theme_void() 

ggsave("analyses/difference.png", p, device = ragg::agg_png)



