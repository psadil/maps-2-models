library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)

d_var <- function(N, d){
  h <- correct_d(N)
  ((N - 1) * (1 + N * d^2) / (N * (N - 3)) - d^2 / h^2) * h^2
}

d_sd <- function(N,d){
  sqrt(d_var(N,d))
}

correct_d <- function(N){
  # https://doi.org/10.1101/865881
  # Han Bossier1∗, Thomas E. Nichols2 & Beatrijs Moerkerke1
  exp((lgamma((N-1)/2)) - log(sqrt((N-1)/2)) - lgamma((N-2)/2))
}


d0 <- tibble(truth = rnorm(1000,0)) |>
  crossing(i = seq_len(100), N=c(10)) |>
  mutate(tstat = rt(n(), N, truth*sqrt(N)) / sqrt(N)) 


quants <- d0 |>
  distinct(truth) |>
  crossing(i = seq_len(1000), N=10) |>
  mutate(tstat = rt(n(), N, truth*sqrt(N)) / sqrt(N)) |>
  group_by(truth, N) |>
  summarise(
    lower = quantile(tstat, 0.05),
    upper = quantile(tstat, 0.95),
    .groups = "drop") |>
  pivot_longer(c(lower, upper))


d0 |>
  group_by(truth, N) |>
  summarise(
    estimate = mean(tstat),
    .groups = "drop") |>
  mutate(estimate = estimate * correct_d(N)) |>
  ggplot(aes(x=truth, y=estimate)) +
  geom_abline() +
  geom_point(alpha=0.1) +
  coord_fixed() +
  geom_line(
    data = quants,
    aes(x=truth, y=value, group=name),
    linetype="dashed"
  )

  



d0 |>
  mutate(estimate = tstat * correct_d(N)) |>
  ggplot(aes(x=truth, y=tstat)) +
  geom_abline() +
  geom_point(alpha=0.1)
  ggpointdensity::geom_pointdensity(show.legend=FALSE, adjust=0.05) +
  scale_color_viridis_c(option="turbo")


d0 |>
  mutate(g = cut(truth, breaks=11)) |>
  group_by(g, N, i) |>
  summarise(
    truth = mean(truth),
    estimate = mean(tstat),
    .groups = "drop") |>
  mutate(estimate = estimate * correct_d(N)) |>
  ggplot(aes(x=truth, y=estimate)) +
  geom_boxplot() +
  geom_abline()


d0 <- tibble(truth = rnorm(1000)) |>
  crossing(i = seq_len(100), N=seq_len(10)) |>
  mutate(x = rnorm(n(), truth)) |>
  group_by(i, truth) |>
  summarise(
    x_bar = mean(x),
    S = sd(x),
    .groups = "drop") 

d0 |>
  mutate(d_hat = (x_bar / S) * correct_d(10)) |>
  group_by(truth) |>
  summarise(estimate = mean(d_hat)) |>
  ggplot(aes(x=truth, y=estimate)) +
  geom_abline() +
  ggpointdensity::geom_pointdensity() +
  scale_color_viridis_c(option="turbo")  +
  coord_fixed()



d0 <- tibble(truth = rnorm(1000)) |>
  crossing(i = seq_len(100), N=seq_len(10)) |>
  mutate(x = rnorm(n(), truth)) |>
  group_by(i, truth) |>
  summarise(
    x_bar = mean(x),
    S = sd(x),
    .groups = "drop") 
d0 |>
  mutate(
    d_hat = (x_bar / S) * correct_d(10)) |>
  group_by(truth) |>
  summarise(
    estimate = mean(d_hat)) |>
  ggplot(aes(x=truth, y=estimate)) +
  geom_abline() +
  ggpointdensity::geom_pointdensity() +
  scale_color_viridis_c(option="turbo")  +
  coord_fixed()



d0 <- tibble(truth = rnorm(1000)) |>
  mutate(observed = rnorm(n(),truth,0.5)) |>
  crossing(i = seq_len(100), N=seq_len(10)) |>
  mutate(x = rnorm(n(), observed)) |>
  group_by(i, truth, observed) |>
  summarise(
    x_bar = mean(x),
    S = sd(x),
    .groups = "drop") 
d0 |>
  mutate(d_hat = (x_bar / S) * correct_d(10)) |>
  group_by(truth) |>
  summarise(estimate = mean(d_hat)) |>
  ggplot(aes(x=truth, y=estimate)) +
  geom_abline() +
  ggpointdensity::geom_pointdensity() +
  scale_color_viridis_c(option="turbo")  +
  coord_fixed()


