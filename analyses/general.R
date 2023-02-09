
# distances

n_pop <- 8526

targets::tar_read(space, store=here::here("_tfce")) |> 
  na.omit() |>
  filter(corrp_thresh==0.95) |>
  whoppeR::WISEsummary(
    dependentvars = "d",
    betweenvars = c("n_sub")
  )

targets::tar_read(space, store=here::here("_tfce")) |> 
  na.omit() |>
  filter(corrp_thresh%in%c(0.95, 0.1)) |> 
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


targets::tar_read(space, store=here::here("_hcp")) |> 
  select(n_sub, Task, corrp_thresh, iter, x, y, z, d) |>
  filter(corrp_thresh%in%c(0.95, 0.1)) |> 
  whoppeR::WISEsummary(
    dependentvars = "d",
    betweenvars = c("n_sub", "corrp_thresh", "Task"),
    na.rm = TRUE
  )
