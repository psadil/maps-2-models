
library(tidyverse)
library(patchwork)

e <- new.env()
tar_load(starts_with("comparison_"), envir = e)

tar_load(comparison)

# comparison |>
p1 <- bind_rows(as.list(e)) |>
  filter(cluster > 0) |>
  mutate(
    pop_avg = pop_avg / sqrt(n_study),
    `N Study` = factor(n_study)) |>
  ggplot(aes(x=pop_avg, y = avg)) +
  geom_point(aes(color = `N Study`)) +
  geom_abline(slope=1, intercept = 0) +
  coord_fixed(xlim = c(0, 5), ylim = c(0,5)) +
  facet_wrap(~n_sub, labeller = "label_both") +
  xlab("average Cohen's d (gold-standard)") +
  ylab("average Z-statistic (ALE)")


p2 <- bind_rows(as.list(e)) |>
  filter(!is.na(dice), dice < .9) |>
  ggplot(aes(x=n_study, y = dice)) +
  facet_wrap(~n_sub, labeller = "label_both") +
  geom_point(alpha = 0.5) +
  scale_y_continuous(
    name = "Dice Coefficient",
    limits = c(0, 1)) +
  xlab("N Studies") 


out <- p1 / p2

ggsave("comparisons.png", out, device = ragg::agg_png)
