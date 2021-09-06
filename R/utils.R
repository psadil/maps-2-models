sample_ukb <- function(d, n, i){
  d |>
    dplyr::slice_sample(n = {{n}}, replace = FALSE) |>
    dplyr::mutate(
      i = {{i}},
      n = {{n}})
}

correlate <- function(d, y){
  d |> 
    dplyr::group_by(i, n) |>
    dplyr::summarise(
      across(tidyselect::matches("[[:digit:]]+"), ~cor({{y}}, .x, method = "spearman")),
      .groups = "drop") |>
    dplyr::select(-{{y}})
}

include_keys <- function(d, key){
  
  d |> 
    tidyr::pivot_longer(
      tidyselect::matches("[[:digit:]]"),
      names_to = "FieldID",
      values_to = "correlation",
      names_transform = list(FieldID = as.integer)) |>
    dplyr::left_join(key, by = "FieldID")
}

load_lit <- function(fname){
  readxl::read_excel(fname) |>
    dplyr::filter(stringr::str_detect(type, "correlation")) |>
    dplyr::select(-notes, -source, -location) |>
    na.omit()
}

plot_rough <- function(d, lit){
  
  p <- d |> 
    ggplot2::ggplot() + 
    ggplot2::facet_grid(
      label~n, labeller = ggplot2::label_both) + 
    ggplot2::geom_histogram(
      ggplot2::aes(x=correlation),
      bins = 50) 
  
  return(p)
}

plot_rough2 <- function(d, lit){
  
  p <- d |> 
    dplyr::mutate(
      N = n,
      Field = stringr::str_remove(Field, " \\(from T1 and T2_FLAIR images\\)")) |>
    ggplot2::ggplot() + 
    ggplot2::facet_wrap(
      ~Field,
      labeller = ggplot2::label_wrap_gen(width = 40)) + 
    ggplot2::geom_point(
      data = lit,
      ggplot2::aes(
        x = value, 
        y = n),
      fill = "green",
      color = "black",
      alpha = 0.3,
      shape = 21) +
    # ggplot2::geom_segment(
    #   data = lit, 
    #   ggplot2::aes(
    #     x = value, 
    #     xend = value, 
    #     y = n, 
    #     yend = 10*n), 
    #   # size=1, 
    #   color = "green") +
    ggridges::stat_density_ridges(
      ggplot2::aes(
        x = correlation, 
        group = N, 
        y = N, 
        fill = factor(stat(quantile))),
      geom = "density_ridges_gradient",
      bandwidth = 0.03,
      calc_ecdf = TRUE, 
      quantiles = c(0.025, 0.975),
      size = 0.01,
      rel_min_height = 0.01,
      scale = 4) + 
    ggplot2::scale_fill_manual(
      name = "Probability",
      values = c("#0000FFA0", "#A0A0A0A0","#FF0000A0"),
      labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")) +
    ggplot2::scale_y_log10() + 
    ggplot2::ggtitle("correlation (Spearman's) with fluid intelligence") +
    ggplot2::theme(legend.position = "bottom") 
  return(p)
}
