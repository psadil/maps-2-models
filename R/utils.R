sample_ukb <- function(d, n, i){
  d |>
    dplyr::slice_sample(n = {{n}}, replace = FALSE) |>
    dplyr::mutate(
      i = {{i}},
      n = {{n}})
}

split_and_correlate <- function(d, n, i, y="fluid"){
  sample_ukb(d=d, n=n, i=i) |> 
    dplyr::group_by(i, n) |>
    dplyr::summarise(
      dplyr::across(
        tidyselect::matches("[[:digit:]]+"), 
        ~cor({{y}}, .x, method = "spearman")),
      .groups = "drop") 
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


dr <- function(r, rho=0.175, n, log = TRUE){
 
  tmp <- purrr::map2_dbl(r, n, ~BAS::hypergeometric2F1(1/2, 1/2, 1/2*(2*.y-1), 1/2*(rho*.x+1), log=TRUE))
  
  if(!log){
    (n-2) * gamma(n-1) * (1 - rho^2)^((n-1)/2) * (1-r^2)^((n-4)/2) / 
      (sqrt(2*pi)*gamma(n-1/2)*(1-rho*r)^(n-3/2)) * exp(tmp)
    # BAS::hypergeometric2F1(1/2, 1/2, 1/2*(2*n-1), 1/2*(rho*r+1), log=FALSE)
  }
  (log(n-2) + lgamma(n-1) + log((1 - rho^2)^((n-1)/2)) + log((1-r^2)^((n-4)/2)) - 
      (log(sqrt(2*pi)) + lgamma(n-1/2) + log((1-rho*r)^(n-3/2))) + tmp) |>
    exp()
  # BAS::hypergeometric2F1(1/2, 1/2, 1/2*(2*n-1), 1/2*(rho*r+1), log=TRUE)
}

get_cor_sampling_distr <- function(n_min=10, n_max=500, rho=0.175){
  
  tidyr::crossing(x = seq(-0.99, .99, length.out=1000), n = seq(n_min, n_max, length.out = 1000)) |>
    dplyr::mutate(y = dr(x, rho=rho, n=n, log=FALSE)) |>
    dplyr::filter(is.finite(y)) |>
    na.omit() |>
    dplyr::group_by(n) |>
    dplyr::arrange(x) |>
    dplyr::mutate(p = cumsum(y)/sum(y)) |>
    dplyr::group_nest() |>
    dplyr::mutate(
      lower = purrr::map_dbl(data, ~max(.x[which.min(abs(0.025 - .x$p)),"x"])),
      upper = purrr::map_dbl(data, ~max(.x[which.min(abs(0.975 - .x$p)),"x"]))) |>
    dplyr::select(-data) |>
    tidyr::pivot_longer(c(lower, upper))
}