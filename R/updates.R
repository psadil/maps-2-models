calc_dice_whole <- function(d, value, lower = 0.0001){
  d |> 
    mutate(
      abs_v = abs({{value}}),
      abs_pop = abs(pop)) |>
    dplyr::group_by(n_sub, iter, n_study) |> 
    dplyr::summarise(
      dice = 2*sum((abs_v > lower) * (abs_pop > lower)) / (sum(abs_v > lower) + sum(abs_pop > lower)),
      .groups = "drop")
}

calc_tpr_whole <- function(d, value, lower = 0.0001){
  d |> 
    mutate(
      abs_v = abs({{value}}),
      abs_pop = abs(pop)) |>
    dplyr::group_by(n_sub, iter, n_study) |> 
    dplyr::summarise(
      dice = sum((abs_v > lower) * (abs_pop > lower)) / sum(abs_pop > lower),
      .groups = "drop")
}

calc_fpr_whole <- function(d, value, lower = 0.1){
  d |> 
    mutate(
      abs_v = abs({{value}}),
      abs_pop = abs(pop)) |>
    dplyr::group_by(n_sub, iter, n_study) |> 
    dplyr::summarise(
      dice = sum((abs_v > lower) * (abs_pop < lower)) / sum(abs_pop < lower),
      .groups = "drop")
}
