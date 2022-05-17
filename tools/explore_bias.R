
library(dplyr)
library(ggplot2)
library(tidyr)

to_tbl0 <- function(value, measure="value"){
  dimnames(value) <- list(
    "x" = seq_len(dim(value)[1]),
    "y" = seq_len(dim(value)[2]),
    "z" = seq_len(dim(value)[3]))
  cubelyr::as.tbl_cube(value, met_name=measure) |> 
    tibble::as_tibble()  
}

to_tbl <- function(file, measure="value", volumes = NULL){
  value <- RNifti::readNifti(file, volumes = volumes)
  to_tbl0(value, measure=measure)
}

get_b <- function(x, y){
  sxx <- sum((x - mean(x))^2)
  syy <- sum((y - mean(y))^2)
  sxy <- sum((x - mean(x)) * (y - mean(y)))
  
  (-1*(sxx - syy) + sqrt((sxx - syy)^2 + 4 * sxy^2))  / (2*sxy)
}

d_var <- function(N, d){
  h <- correct_d(N)
  ((N - 1) * (1 + N * d^2) / (N * (N - 3)) - d^2 / h^2) * h^2
}

correct_d <- function(N){
  # https://doi.org/10.1101/865881
  # Han Bossier1∗, Thomas E. Nichols2 & Beatrijs Moerkerke1
  exp((lgamma((N-1)/2)) - log(sqrt((N-1)/2)) - lgamma((N-2)/2))
}

fix_names <- function(d){
  d |>
    dplyr::mutate(
      tstat = stringr::str_replace(.data$tstat, "/_tstat1.nii.gz", "_glm_sigmasqr.nii.gz")
    )
}

mask <- function(d, mask=MNITemplate::getMNIPath("Brain_Mask", "2mm")){
  m <- to_tbl0(RNifti::readNifti(mask)) |>
    dplyr::filter(value > 0)
  dplyr::semi_join(d, m, by = c("x","y","z"))
}

do_cor0 <- function(studies){
  studies |>
    fix_names() |>
    dplyr::mutate(
      rho = purrr::map(
        .data$tstat, 
        ~to_tbl(.x) |> mask()
      )
    )
}


pop_cope <- to_tbl("data-raw/niis/nsub-8526_iter-0_glm_cope.nii.gz") |>
  mask()
cope500 <- to_tbl("data-raw/niis/nsub-500_iter-1_glm_cope.nii.gz") |>
  mask()

left_join(pop_cope, cope500, by = c("x","y","z")) |>
  slice_sample(n=100000) |>
  ggplot(aes(x=value.x, y=value.y)) +
  geom_point(alpha = 0.02) +
  coord_fixed() +
  geom_abline(slope=1, intercept=0) +
  geom_smooth(method = "lm")

pop_var <- to_tbl("data-raw/niis/nsub-8526_iter-0_glm_sigmasqr.nii.gz") |>
  mask()
var500 <- to_tbl("data-raw/niis/nsub-500_iter-3_glm_sigmasqr.nii.gz") |>
  mask()

left_join(pop_var, var500, by = c("x","y","z")) |>
  slice_sample(n=100000) |>
  mutate(across(starts_with("value"), sqrt)) |>
  ggplot(aes(x=value.x, y=value.y)) +
  geom_point(alpha = 0.02) +
  coord_fixed() +
  geom_abline(slope=1, intercept=0) +
  geom_smooth(method = "lm")

pop_tstat <- to_tbl("data-raw/niis/nsub-8526_iter-0_tstat1.nii.gz") |>
  mutate(value = value / sqrt(8526) * correct_d(8526))
tstat500 <- to_tbl("data-raw/niis/nsub-500_iter-1_tstat1.nii.gz") |>
  mutate(value = value / sqrt(500) * correct_d(500))

left_join(pop_tstat, tstat500, by = c("x","y","z")) |>
  filter(abs(value.x) > 0.01) |>
  slice_sample(n=100000) |>
  ggplot(aes(x=value.x, y=value.y)) +
  geom_point(alpha = 0.02) +
  coord_fixed() +
  geom_abline(slope=1, intercept=0) +
  geom_smooth(method = "lm")


tmp <- do_cor0(filter(tfce, n_sub==10)) |>
  select(iter, n_sub, rho) |> 
  unnest(rho) |>
  right_join(pop_var, by = c("x","y","z"), suffix = c(".study",".pop")) 

tmp |>
  # filter(value.pop > 50^2) |>
  slice_sample(n=100000) |>
  mutate(across(starts_with("value"), sqrt)) |>
  ggplot(aes(x=value.pop, y=value.study)) +
  geom_point(alpha = 0.02) +
  # coord_fixed() +
  geom_abline(slope=1, intercept=0) +
  geom_smooth(method = "lm")

tmp |>
  group_by(x, y, z) |>
  summarise(
    value.pop = mean(value.pop), 
    value.study = mean(value.study), 
    .groups = "drop") |>
  # slice_sample(n=100000) |>
  mutate(across(starts_with("value"), sqrt)) |>
  ggplot(aes(x=value.pop, y=value.study)) +
  geom_point(alpha = 0.02) +
  # coord_fixed() +
  geom_abline(slope=1, intercept=0) +
  geom_smooth(method = "lm")



tmp <- targets::tar_read(tfce_cor) |>
  select(iter, n_sub, d) |>
  unnest(d) 


tmp2 <- tmp |>
  group_by(x, y, z, n_sub) |>
  summarise(
    value = mean(value), 
    .groups = "drop") |>
  right_join(
    to_tbl(tfce_pop$tstat) |> mask() |> mutate(value = value / sqrt(8526) * correct_d(8526)), 
    by = c("x","y","z"), 
    suffix = c(".study",".pop"))

tmp2 |>
  # filter(value.pop > 50^2) |>
  # slice_sample(n=100000) |>
  ggplot(aes(x=value.pop, y=value.study)) +
  geom_point(alpha = 0.02) +
  facet_wrap(~n_sub) +
  # coord_fixed() +
  geom_abline(slope=1, intercept=0) +
  geom_smooth()


tmp <- targets::tar_read(tfce_cor) |>
  select(iter, n_sub, d) |>
  unnest(d) 

targets::tar_load(tfce_pop) 
pop_tstat <- to_tbl(tfce_pop$tstat, measure = "gold") |> 
  mask() |> 
  mutate(
    gold = gold / sqrt(8526) * correct_d(8526),
    group = case_when(
      abs(gold) < .2 ~ "nill",
      between(abs(gold), .2, .5) ~ "small",
      between(abs(gold), .5, .8) ~ "medium",
      abs(gold) > .8 ~ "large",
    ),
    group = factor(group, levels = c("nill","small","medium","large"))) 

tmp2 <- tmp |>
  filter(n_sub==10) |>
  group_by(x, y, z, n_sub) |>
  summarise(
    value = mean(value), 
    .groups = "drop") |>
  left_join(pop_tstat)

tmp2 |>
  # filter(gold > .5) |>
  # slice_sample(n=100000) |>
  ggplot(aes(x=gold, y=value)) +
  facet_wrap(~n_sub) +
  geom_point(alpha = 0.02) +
  geom_abline(slope=1, intercept=0) +
  geom_smooth() 


iters <- tmp |>
  right_join(
    filter(pop_tstat, !forcats::fct_match(group, "nill")) |>
      mutate(gold_cut = cut(gold, breaks=100))) |>
  group_by(gold_cut, n_sub, iter) |>
  summarise(
    s = sd(value),
    N = n(),
    q.05 = quantile(value, 0.05),
    q.50 = median(value),
    q.95 = quantile(value, 0.95),
    gold = mean(gold),
    .groups = "drop") |>
  group_by(gold_cut, n_sub) |>
  summarise(
    s = median(s),
    N = median(N),
    q.05 = median(q.05),
    q.50 = median(q.50),
    q.95 = median(q.95),
    gold = median(gold),
    .groups = "drop")

at <- qs::qread("data-raw/atlas.qs")


filter_space <- function(d, n_peaks, store="_tfce"){
  targets::tar_read(space, store=store) |>
    # filter(between(z, 25, 70)) |>
    # filter(`Cluster Index`==10) |>
    filter(!is.na(x.study)) |>
    group_nest() |>
    mutate(top = map(data, ~distinct(.x, Value) |> slice_max(order_by=Value, n=.env$n_peaks))) |>
    unnest(data) |>
    mutate(keep = map2_lgl(Value, top, ~any(.x %in% .y$Value))) |>
    filter(keep) |>
    mutate(n_sub = factor(n_sub), group = interaction(x,y,z, lex.order = TRUE))
}

filtered_space <- filter_space(space, 9) |>
  distinct(label) |>
  filter(str_detect(label, "Amyg|Occipital Fusiform|Occipital Pole")) |>
  mutate(
    text = case_when(
      label == "Temporal Occipital Fusiform Cortex" ~ "Temp. Occipital Fus.",
      label == "Occipital Fusiform Gyrus" ~ "Occ. Fus. Gyr.",
      str_detect(label, "Amyg")  ~ "Amygdala",
      TRUE ~ "Occipital Pole"
    ))

gold_meds <- pop_tstat |>
  left_join(at) |>
  right_join(filtered_space) |>
  group_by(text) |>
  summarise(gold = median(gold)) |>
  mutate(
    y = c(-.5, -.75, -1, -.25),
    hjust = c("left", "right", "right", "left"))

# not using, because not clear why 90% bounds so poor for low n_sub
# theoretical <- iters |>
#   distinct(gold, n_sub) |>
#   mutate(
#     q.05 = qt(.05, n_sub-1, gold*sqrt(n_sub))/sqrt(n_sub),
#     q.50 = qt(.50, n_sub-1, gold*sqrt(n_sub))/sqrt(n_sub),
#     q.95 = qt(.95, n_sub-1, gold*sqrt(n_sub))/sqrt(n_sub),
#     src = "Assumption Based")

sigs <- distinct(iters, n_sub) |>
  crossing(p = c(0.01, 0.001)) |>
  mutate(
    lower = qt(p, n_sub-1) / sqrt(n_sub),
    upper = qt(p, n_sub-1, lower.tail = FALSE) / sqrt(n_sub))

iters |>
  mutate(src = "empirical") |>
  # bind_rows(theoretical) |>
  select(-s, -average) |>
  pivot_longer(cols = c(q.05, q.95, q.50), names_to = "stat") |>
  mutate(quantile = if_else(stat=="q.50", "median", "90%")) |>
  ggplot(aes(x=gold)) +
  facet_wrap(~n_sub) +
  geom_vline(
    data = gold_meds,
    aes(xintercept = gold),
    alpha = 0.25) +
  geom_text(
    data = gold_meds,
    aes(x = gold, y=y, label=text, hjust=hjust)
  ) +
  geom_hline(
    data = sigs,
    aes(yintercept = upper),
    alpha = 0.25) +
  geom_text(
    data = mutate(
      sigs, 
      upper = if_else(p==0.01, upper-.1, upper+.1),
      p = glue::glue("p = {p}")),
    aes(y=upper, label = p),
    x = -0.5) +
  geom_abline(slope=1, intercept=0, alpha=0.5) +
  geom_line(aes(y=value, linetype = quantile, group=stat)) +
  scale_linetype_manual(
    name = "Study Quantile",
    breaks = c("median", "90%"),
    values = c("median"="dashed", "90%"="dotted")) +
  xlab("Gold Standard d") +
  ylab("Study d")
  



# purrr::walk(train[1:100], ~fs::file_copy(.x, "tools/tmp"))
manual <- tibble(cope = fs::dir_ls("tools/tmp")) |>
  slice_head(n=10) |>
  mutate(data = map(cope, to_tbl)) |>
  unnest(data) |>
  mask() |>
  group_by(x, y, z) |>
  summarise(
    s = var(value),
    mu = mean(value),
    N = n(),
    .groups = "drop")

fsld <- to_tbl("tools/rando__glm_sigmasqr.nii.gz") |>
  mask() 

left_join(manual, fsld) |>
  mutate(di = s - value) |>
  summarise(
    m = mean(di),
    s = sd(di)
  )
