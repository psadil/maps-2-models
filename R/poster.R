
filter_space <- function(d, n_peaks){
  d |>
    # filter(between(z, 25, 70)) |>
    # filter(`Cluster Index`==10) |>
    dplyr::filter(!is.na(x.study)) |>
    dplyr::group_nest() |>
    dplyr::mutate(
      # first take the top from each region
      top = purrr::map(
        data, 
        ~dplyr::distinct(.x, Value, x, y, z, label, hemi) |> 
          dplyr::group_by(label, hemi) |>
          dplyr::slice_max(order_by=Value, n=1)),
      data = purrr::map2(
        data, top, 
        ~dplyr::semi_join(.x, .y, by = c("x","y","z","Value"))),
      # then the top n_peaks
      top = purrr::map(
        data, 
        ~dplyr::distinct(.x, Value, x, y, z) |> 
          dplyr::slice_max(order_by=Value, n=.env$n_peaks)),
      data = purrr::map2(
        data, top, 
        ~dplyr::semi_join(.x, .y, by = c("x","y","z","Value")))) |>
    dplyr::select(-top) |>
    tidyr::unnest(data) |>
    dplyr::mutate(
      nsub = factor(n_sub), 
      group = interaction(x,y,z, lex.order = TRUE))
}


make_big_plot <- function(filename, big, iters, space,pop_d, at){
  
  sigs <- dplyr::distinct(iters, n_sub) |>
    tidyr::crossing(p = c(0.01, 0.001)) |>
    dplyr::mutate(
      lower = qt(p, n_sub-1) / sqrt(n_sub),
      upper = qt(p, n_sub-1, lower.tail = FALSE) / sqrt(n_sub)) |>
    dplyr::rename(`N Participants Simulated` = n_sub)
  
  filtered_space <- filter_space(space, 9) |>
    dplyr::distinct(label) |>
    dplyr::filter(stringr::str_detect(label, "Amyg|Occipital Fusiform Gyrus|Occipital Pole")) |>
    dplyr::mutate(
      text = dplyr::case_when(
        label == "Occipital Fusiform Gyrus" ~ "Occ. Fus. Gyr.",
        stringr::str_detect(label, "Amyg")  ~ "Amygdala",
        TRUE ~ "Occipital Pole"
      ))
  
  gold_meds <- pop_d |>
    dplyr::left_join(at) |>
    dplyr::right_join(filtered_space) |>
    dplyr::group_by(text) |>
    dplyr::summarise(gold = median(gold)) 
  
  # targets::tar_load(tfce_cor, store = params$tfce_store)
  # a0 <- tfce_cor$d[[310]] |>
  #   slice_sample(prop=.01) |>
  #   left_join(gold0) |>
  #   ggplot(aes(x=gold, y=value)) +
  #   ggpointdensity::geom_pointdensity(show.legend = FALSE)  +
  #   scale_color_viridis_c(option="turbo", name=NULL) +
  #   # coord_fixed() +
  #   geom_abline() +
  #   geom_vline(
  #     xintercept=seq(min(gold0$gold), max(gold0$gold), length.out=11),
  #     alpha = 0.5) +
  #   scale_x_continuous(
  #     name = "Gold standard d",
  #     breaks = c(-1, 0, 1),
  #     labels = c(-1, 0, 1),
  #     limits = c(-1, 2)) +
  #   scale_y_continuous(
  #     name = bquote(paste(Study[i], " g")),
  #     breaks = c(-1, 0, 1),
  #     labels = c(-1, 0, 1),
  #     limits = c(-1.5, 2)) +
  #   theme(
  #     axis.ticks = element_blank(),
  #     axis.text.x = element_blank(),
  #     axis.text.y = element_blank(),
  #     axis.title.y = element_text(size = 12, margin = margin(r = -20)),
  #     axis.title.x = element_text(size = 12, margin = margin(t= -20)),
  #     plot.background = element_blank(),
  #     panel.grid = element_blank())
  # 
  # a1 <- a0 + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  draw_line <- function(N,q, lower=min(big$gold), upper=max(big$gold)){
    d <- tibble::tibble(`N Participants Simulated`=N, gold=seq(lower, upper, length.out=10))
    ggplot2::geom_function(
      data = d,
      fun = ~qt(q, N-1, .x*sqrt(N))/sqrt(N), linetype="dashed")
  }
  
  p <- big |>
    # dplyr::filter(n_sub==10) |>
    dplyr::rename(`N Participants Simulated` = n_sub) |>
    dplyr::group_by(`N Participants Simulated`) |>
    dplyr::slice_sample(prop=0.01) |>
    ggplot2::ggplot(ggplot2::aes(x=gold)) +
    ggplot2::facet_wrap(~`N Participants Simulated`, labeller = ggplot2::label_both) +
    # geom_point(aes(y=value), alpha=0.1) +
    ggplot2::geom_density_2d_filled(
      ggplot2::aes(y=value), 
      breaks = c(.01, seq(.1, 1, by=.1)),
      contour_var = "ndensity") +
    # ggpointdensity::geom_pointdensity(
    #   ggplot2::aes(y=value), 
    #   show.legend = FALSE, 
    #   adjust = 0.01) +
    ggplot2::scale_fill_viridis_d(option="turbo", name="Normalized\nDensity") +
    # ggnewscale::new_scale_fill() +
    draw_line(10, 0.95) +
    draw_line(10, 0.05) +
    draw_line(20, 0.95) +
    draw_line(20, 0.05) +
    draw_line(50, 0.95) +
    draw_line(50, 0.05) +
    draw_line(100, 0.95) +
    draw_line(100, 0.05) +
    draw_line(200, 0.95) +
    draw_line(200, 0.05) +
    draw_line(500, 0.95) +
    draw_line(500, 0.05) +
    ggplot2::geom_vline(
      data = gold_meds,
      ggplot2::aes(xintercept = gold, color=text),
      size = 2) +
    ggplot2::scale_color_viridis_d(option="turbo", name=NULL, begin=0.1) +
    ggplot2::geom_hline(
      data = sigs,
      ggplot2::aes(yintercept = upper),
      alpha = 0.25,
      size = 1) +
    ggplot2::geom_text(
      data = dplyr::mutate(
        sigs,
        upper = dplyr::if_else(p==0.01, upper-.2, upper+.2),
        p = glue::glue("p = {p}")),
      ggplot2::aes(y=upper, label = p), x = -1, size = 5, hjust="left") +
    ggplot2::geom_abline(color="darkgoldenrod1") +
    ggplot2::scale_x_continuous(
      name = "Gold standard d",
      breaks = c(-1, 0, 1),
      labels = c(-1, 0, 1),
      limits = c(-1, 2)) +
    ggplot2::scale_y_continuous(
      name = "Study g",
      breaks = c(-2, 0, 2),
      labels = c(-2, 0, 2),
      limits = c(-2, 2)) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.justification = "left",
      legend.key.size = ggplot2::unit(0.5, "in"),
      legend.spacing.x = ggplot2::unit(0.3, "in")
    )
  
  
  
  # b <- iters |>
  #   rename(`N Participants Simulated` = n_sub) |>
  #   ggplot(aes(x=gold)) +
  #   facet_wrap(~`N Participants Simulated`, labeller = label_both) +
  #   geom_vline(
  #     data = gold_meds,
  #     aes(xintercept = gold, color=text),
  #     size = 2) +
  #   geom_hline(
  #     data = sigs,
  #     aes(yintercept = upper),
  #     alpha = 0.25,
  #     size = 1) +
  #   geom_text(
  #     data = mutate(
  #       sigs,
  #       upper = if_else(p==0.01, upper-.2, upper+.2),
  #       p = glue::glue("p = {p}")),
  #     aes(y=upper, label = p),
  #     x = -0.5,
  #     size = 5) +
  #   geom_abline(alpha=0.5, size=1,color="darkgoldenrod1") +
  #   geom_boxplot(aes(color=gold, y=average)) +
  #    scico::scale_color_scico(palette = "vik") +
  #   # geom_line(aes(y=value, linetype = quantile, group=stat),  size=1.5) +
  #   scale_x_continuous(
  #     name = "Gold standard d",
  #     breaks = c(-1, 0, 1),
  #     labels = c(-1, 0, 1),
  #     limits = c(-1, 2)) +
  #   scale_y_continuous(
  #     name = "avg(Study g)",
  #     breaks = c(-1, 0, 1),
  #     labels = c(-1, 0, 1),
  #     limits = c(-1.5, 2)) +
  #   # scale_linetype_manual(
  #   #   name = "Study\nStatistic",
  #   #   breaks = c("average", "90%"),
  #   #   values = c("average"="dotted", "90%"="dashed")) +
  #   # scale_color_viridis_d(option="turbo", begin = .1, name=NULL) +
  #   theme(
  #     legend.position = "bottom",
  #     legend.justification = "left",
  #     legend.key.size = unit(0.5, "in"),
  #     legend.spacing.x = unit(0.3, "in")
  #     # legend.key = element_rect(color = NA, fill = NA)
  #   )
  
  # a1 + a1 + a1 + (a0+ggtitle("Equal Bin Width")) + b +
  #   patchwork::plot_layout(
  #     design=c(
  #       area(t = 4, l = 4, b = 8, r = 8),
  #       area(t = 3, l = 3, b = 7, r = 7),
  #       area(t = 2, l = 2, b = 6, r = 6),
  #       area(t = 1, l = 1, b = 5, r = 5),
  #       area(t = 1, l = 9, b = 8, r = 50)
  #     )) 
  
  # ggplot2::ggsave(
  #   filename = filename, 
  #   p, 
  #   device = ragg::agg_png(
  #     res=300,
  #     units = "in",
  #     height=7, 
  #     width=18
  #   ))  
  # 
  p
}


# https://themockup.blog/posts/2020-08-28-heatmaps-in-ggplot2/#get-the-raw-density-estimates
get_density <- function(x, y, ...) {
  # density_out <- MASS::kde2d(x, y, ...)
  h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))
  density_out <- KernSmooth::bkde2D(cbind(x, y), bandwidth=h, gridsize=c(100, 100), ...)
  # int_x <- findInterval(x, density_out$x)
  # int_y <- findInterval(y, density_out$y)
  
  int_x <- findInterval(x, density_out$x1)
  int_y <- findInterval(y, density_out$x2)
  comb_int <- cbind(int_x, int_y)
  # return(density_out$z[comb_int])
  return(density_out$fhat[comb_int])
}


make_big <- function(tfce_cor, pop_d, prop=0.01){
  
  if(!all.equal(prop, 1)){
    tmp <- pop_d |>
      dplyr::slice_sample(prop = prop)
  }else{
    tmp <- pop_d
  }
  
  tfce_cor |>
    dplyr::select(iter, n_sub, d) |>
    tidyr::unnest(d) |>
    dplyr::right_join(tmp, by = c("x","y","z")) |>
    dplyr::select(x,y,z,value,gold,iter,n_sub) |>
    dplyr::group_nest(n_sub) |>
    dplyr::mutate(
      data = purrr::map(
        data,
        ~dplyr::mutate(.x, dens = get_density(gold, value)))) |>
    tidyr::unnest(data)
}

get_gold_bvf_cor <- function(ukb, test){
  ukb_test <- dplyr::semi_join(
    ukb, 
    tibble::tibble(
      f.eid = stringr::str_extract(test, "[[:digit:]]{7}") |> as.integer()), 
    by = "f.eid")
  cor(ukb_test$`25010`, ukb_test$fluid, method = "spearman")
}
