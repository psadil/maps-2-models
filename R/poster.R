
filter_space <- function(d, n_peaks){
  d |>
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


