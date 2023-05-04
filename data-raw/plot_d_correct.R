library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)


d_var <- function(N, d){
  h <- correct_d(N)
  ((N - 1) * (1 + N * d^2) / (N * (N - 3)) - d^2 / h^2) * h^2
}

correct_d <- function(N){
  # https://doi.org/10.1101/865881
  # Han Bossier1∗, Thomas E. Nichols2 & Beatrijs Moerkerke1
  # https://journals.sagepub.com/doi/pdf/10.3102/10769986006002107?casa_token=d0UqQNTgaDMAAAAA:qN_ot_q6J-K6D6iRgId5qy1Edqldg8xcRZ54k_dhEfiYoywEU6dczwwb_GXkjhXhoWa7jetEfsUQqw
  exp((lgamma((N-1)/2)) - log(sqrt((N-1)/2)) - lgamma((N-2)/2))
}


d <- crossing(N = c(10, 100), i = seq_len(50), d = 2) |>
  mutate(
    copes = map2(N,d, ~rnorm(.x, .y)),
    m = map_dbl(copes, mean),
    s = map_dbl(copes, sd),
    cd = m / s,
    hg = cd * correct_d(N)) |>
  select(-copes) |>
  pivot_longer(cols = c(m, s, cd, hg))

d |>
  ggplot(aes(x = value)) +
  facet_grid(name~N) +
  geom_histogram(bins=100)

d |>
  group_by(name, N) |>
  summarise(
    m = mean(value),
    .groups = "drop")


short <- crossing(N = c(10, 20), d = rnorm(10000, 2)) |>
  mutate(
    cd = rt(n(), N, d*sqrt(N)) / sqrt(N),
    hg = cd * correct_d(N)
  )

short |>
  ggplot(aes(x=d, y=hg)) +
  facet_wrap(~N) +
  coord_fixed() +
  geom_point(alpha = 0.1) +
  geom_abline(intercept = 0, slope=1) +
  geom_smooth(method = "lm")
