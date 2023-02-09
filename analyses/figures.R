library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)
source(here::here("R","utils.R"))
source(here::here("R","tfce.R"))
source(here::here("R","updates.R"))
source(here::here("R","poster.R"))

yeo_dist <- function(filename=here::here("analyses", "agv-dist-yeo-ukb.png")){
  n_pop <- targets::tar_read(tfce_pop, store=here::here("_tfce"))$n_sub
  load_space <- function(store=here::here("_tfce")){
    targets::tar_read(space, store = store) |>
      filter(!is.na(label)) |>
      mutate(corrp_thresh = factor(corrp_thresh)) |>
      filter_space(n_peaks=100) |>
      rename(Network = `Network Name`) |>
      mutate(
        Value = Value / sqrt(n_pop),
        label = if_else(
          is.na(`Network`), 
          label, 
          str_extract(label, "(?<=_)([[:alnum:]]+)_[[:digit:]]+$")),
        `Network` = if_else(
          is.na(`Network`), 
          "subcortical", 
          `Network`
        ),
        `Network` = factor(
          `Network`, 
          levels = c(
            "limbic A",
            "default A", "default B", "default C",
            "control A", "control B",
            "salience / ventral attention A",
            "salience / ventral attention B",
            "dorsal attention A","dorsal attention B",
            "temporal parietal",
            "somatomotor A", "somatomotor B", "central visual", "peripheral visual",
            "subcortical"
          ),
          labels = c(
            "lmbc A",
            "def A", "def B", "def C",
            "cntrl A", "cntrl B",
            "sal A", "sal B",
            "dors atn A", "dors atn B",
            "temp par",
            "sm A", "sm B", "cnt vis", "periph vis",
            "sub cort"
          ),
          ordered = TRUE)) 
  }
  
  space5 <- load_space(here::here("_tfce")) |>
    mutate(contrast = "Face > Shape")
  space1 <- load_space(here::here("_tfce1"))  |>
    mutate(contrast = "Shape > 0") 
  
  space <- bind_rows(space5, space1) |>
    filter(corrp_thresh == "0.95") |>
    mutate(contrast = factor(contrast))
  
  p <- space  |>
    rename(`N Sub` = n_sub) |>
    ggplot(aes(y=`Network`, x=d)) +
    facet_grid(contrast~`N Sub`, labeller = label_both, shrink = TRUE, scales = "free_y") +
    geom_boxplot(outlier.alpha = 0.25) +
    scale_x_continuous(
      name = "avg dist(Gold Peak, Study Peak) (mm)",
      breaks = c(0, 100),
      labels = c(0, 100)) +
    theme_gray(base_size = 10)
  
  ggsave(
    filename = filename,
    plot = p, 
    device = ragg::agg_png,
    height = 4,
    width=5.5)
}



fluid <- function(filename = here::here("analyses", "fluid.png")){
  targets::tar_load(ukb_cor, store = "_tfce")
  brainvolfluid <- targets::tar_read(brainvolfluid, store = "_tfce")
  
  lit <- targets::tar_read(lit, store = "_tfce") |>
    filter(n>4)
  
  targets::tar_load(sampling_distr, store = "_tfce")
  
  s <- 2.5
  p <- brainvolfluid |> 
    mutate(N = factor(n)) |>
    ggplot2::ggplot(aes(y=n, x = correlation)) + 
    geom_vline(xintercept = ukb_cor, color = "darkgoldenrod1") +
    geom_boxplot(
      aes(group=n),
      outlier.shape = NA) +
    geom_jitter(width = 0, alpha = 0.05, size=s) +
    ggplot2::geom_point(
      data = lit,
      ggplot2::aes(x = value, y = n),
      fill = "green",
      color = "black",
      size=s,
      # alpha = 0.3,
      shape = 21) +
    geom_line(
      data = sampling_distr,
      aes(x=value, y=n, group=name)) +
    xlab("cor(Brain Volume, Fluid Intelligence)") + 
    ylab("N Participants in Simulation,\nPublished Study") +
    ggplot2::scale_y_log10(
      labels = c(1, 10, 100, 1000),
      breaks = c(1, 10, 100, 1000),
      limits = c(5, 1000)) +
    annotate(geom="point", fill="green", color="black", shape=21, x=-.85, y=800, size=s) +
    annotate(geom = "text", x=-.8, y=800, label = "Published Correlation", hjust = "left", size=s) +
    annotate(geom="point", x=-.85, y=500, alpha=0.5, size=s) +
    annotate(geom="text", x=-.8, y=500, label = "Simulated Study Correlation", hjust = "left", size=s) +
    annotate(geom="segment", x=-.85, xend=-.85, yend=250, y=350, color="darkgoldenrod1", size=s) +
    annotate(geom="text", x=-.8, y=300, label = "Gold Standard Correlation", hjust = "left", size=s) +
    annotate(geom="segment", x=-.85, xend=-.9, yend=100, y=200) +
    annotate(geom="text", x=-.8, y=140, label = "Theoretical Sampling\nDistribution", hjust = "left", size=s) +
    theme_gray(base_size = 9)
  
  ggsave(
    filename=filename,
    plot = p,
    device = ragg::agg_png,
    width = 4,
    height = 3
  )
}

voxelwise <- function(filename = here::here("analyses", "voxelwise.png")){
  targets::tar_load(iters, store = "_tfce")
  targets::tar_load(pop_d, store = "_tfce")
  
  sigs <- iters |>
    dplyr::distinct(n_sub) |>
    tidyr::crossing(p = c(0.01, 0.001)) |>
    dplyr::mutate(
      lower = qt(p, n_sub-1) / sqrt(n_sub),
      upper = qt(p, n_sub-1, lower.tail = FALSE) / sqrt(n_sub)) |>
    dplyr::rename(`N Participants Simulated` = n_sub)
  
  filtered_space <- filter_space(targets::tar_read(space, store="_tfce"), 9) |>
    dplyr::distinct(label) |>
    dplyr::filter(stringr::str_detect(label, "Amyg|Occipital Fusiform Gyrus|Occipital Pole")) |>
    dplyr::mutate(
      text = dplyr::case_when(
        label == "Occipital Fusiform Gyrus" ~ "Occ. Fus. Gyr.",
        stringr::str_detect(label, "Amyg")  ~ "Amygdala",
        TRUE ~ "Occipital Pole"
      ))
  
  gold_meds <- pop_d |>
    dplyr::left_join(targets::tar_read(at, store = "_tfce")) |>
    dplyr::right_join(filtered_space) |>
    dplyr::group_by(text) |>
    dplyr::summarise(gold = median(gold)) 
  
  
  p <- targets::tar_read(big, store = "_tfce") |>
    dplyr::rename(`N Participants Simulated` = n_sub) |>
    ggplot2::ggplot(ggplot2::aes(x=gold)) +
    ggplot2::facet_wrap(~`N Participants Simulated`, labeller = ggplot2::label_both) +
    scattermore::geom_scattermore(aes(y=value, color=dens/max(dens), alpha=dens)) +
    ggplot2::scale_color_viridis_c(
      option="turbo", 
      name="Voxel\nDensity",
      labels = c(0,1),
      breaks = c(0,1),
      limits=c(0,1)) +
    guides(alpha = "none") + 
    # ggnewscale::new_scale_color() +
    # ggplot2::geom_vline(
    #   data = gold_meds,
    #   ggplot2::aes(xintercept = gold, color=text),
    #   size = 2) +
    # ggplot2::scale_color_viridis_d(option="turbo", name="ROI", begin=0.1) +
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
      ggplot2::aes(y=upper, label = p), x = -1, size = 2, hjust="left") +
    ggplot2::geom_abline(color="darkgoldenrod1") +
    ggplot2::scale_x_continuous(
      name = "Gold Standard Effect Size (d)",
      breaks = c(-1, 0, 1),
      labels = c(-1, 0, 1),
      limits = c(-1, 2)) +
    ggplot2::scale_y_continuous(
      name = "Study Effect Size (g)",
      breaks = c(-2, 0, 2),
      labels = c(-2, 0, 2),
      limits = c(-2, 2)) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.justification = "left",
      legend.key.size = ggplot2::unit(0.2, "in"),
      legend.spacing.x = ggplot2::unit(0.3, "in"),
      legend.text = element_text(size = 18),
      legend.background = element_rect(fill = "transparent"),
    )
  
  ggsave(filename=filename, plot = p, device = ragg::agg_png, width=17.8, units = "cm", height = 4*2.54)
}


gold_pk <- function(filename = here::here("analyses", "gold_peaks.png")){
  
  space <- targets::tar_read(space, store="_tfce") |>
    filter(str_detect(label, "Stem", negate = TRUE)) |>
    filter_space(n_peaks=Inf) |>
    filter(!is.na(`Full Component Name`)) |>
    group_by(`Full Component Name`) |>
    slice_max(order_by=Value, n=1) |>
    ungroup()
  
  zs <- distinct(space, z)$z
  
  targets::tar_load(tfce_pop, store="_tfce")
  
  p <- to_tbl(tfce_pop$tstat[[1]]) |>
    mutate(value = value / sqrt(tfce_pop$n_sub[[1]])) |>
    filter(z %in% zs, abs(value)>0.01) |>
    mutate(a=abs(value)) |>
    ggplot(aes(x=x, y=y)) +
    coord_fixed(clip = "off") +
    facet_wrap(~z, labeller=label_both, ncol=6) +
    geom_raster(
      aes(fill=value),
      data = to_tbl(MNITemplate::getMNIPath(what="Brain", res = "2mm")) |> 
        filter(z%in%zs), show.legend = FALSE) +
    scale_fill_distiller(
      type = "seq",
      direction = -1,
      palette = "Greys") +
    ggnewscale::new_scale_fill() +
    geom_raster(aes(alpha=a, fill=value)) +
    scico::scale_fill_scico(
      palette = "vik",
      limits=c(-2, 2),
      breaks = c(-1.5, 1.5),
      labels = c(-1.5, 1.5)) +
    scale_alpha_continuous(guide="none", range = c(0.5, 1)) +
    geom_point(
      shape = 4,
      color="red", stroke=0.5,
      data=distinct(space, x,y,z), 
      size=2,
      alpha=0.5) +
    # ggrepel::geom_label_repel(
    #   aes(label=label), alpha=0.75,
    #   data=distinct(space, x,y,z, label)
    # ) +
    labs(fill = "Face > Shape\n(Cohen's d)") +
    theme_void(base_size = 9) +
    labs(caption = "<span style='color:red'>X</span> Local Peak") +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
      # legend.spacing.x= unit(1.0, "in"),
      # legend.key.width = unit(0.5, 'in'),
      plot.caption = ggtext::element_markdown(lineheight = 1.2)
    )  
  
  ggsave(filename=filename, plot=p, device=ragg::agg_png, width=6, height = 4)
}


bysize <- function(filename=here::here("analyses", "bysize.png")){
  
  p <- targets::tar_read(space, store=here::here("_tfce")) |> 
    rename(Network = `Network Name`, FWE = corrp_thresh, `N Sub`=n_sub) |>
    filter(!is.na(label)) |>
    mutate(Value = Value / sqrt(n_pop)) |>
    group_by(Value, `N Sub`, FWE, Network, label, hemi) |>
    summarise(d = mean(d, na.rm=TRUE), .groups = "drop") |>
    ggplot(aes(y=d, x=Value)) +
    # geom_line(aes(group=Value), alpha=0.5) +
    # geom_vline(aes(xintercept=upper, linetype=p), data=sigs) +
    # geom_point(aes(color=corrp_thresh), alpha=0.5) +
    geom_point(alpha=0.5) +
    facet_grid(FWE~`N Sub`, labeller = label_both) +
    scale_x_continuous(
      "Gold Standard Peak Cohen's d",
      breaks = c(0, 1),
      labels = c(0, 1)) +
    ylab("avg dist(Gold Peak, Study Peak) (mm)") +
    theme_gray(base_size=9)
  
  ggsave(filename, plot=p, device=ragg::agg_png, height = 3, width=5.5)
}


# HCP

avg_within_network <- function(){
  
  store<- "_hcp"
  
  load_space <- function(){
    targets::tar_read(space, store = store) |>
      filter(!is.na(label)) |>
      mutate(corrp_thresh = factor(corrp_thresh)) |>
      rename(Network = `Network Name`) |>
      mutate(
        Value = Value / sqrt(n_pop),
        label = if_else(
          is.na(`Network`), 
          label, 
          str_extract(label, "(?<=_)([[:alnum:]]+)_[[:digit:]]+$")),
        `Network` = if_else(
          is.na(`Network`), 
          "subcortical", 
          `Network`
        ),
        `Network` = factor(
          `Network`, 
          levels = c(
            "limbic",
            "default",
            "control",
            "salience / ventral attention",
            "dorsal attention",
            "somatomotor", "visual",
            "subcortical"
          ),
          ordered = TRUE)
      ) 
  }
  space <- load_space() |>
    left_join(targets::tar_read(contrasts, store=store),by = c("Task", "CopeNumber")) |>
    rename(contrast=ContrastName, FWE = corrp_thresh) 
  
  p <- space |>
    filter(fct_match(FWE, "0.01")) |>
    group_by(Task, contrast, n_sub, iter, Network) |>
    summarise(d = mean(d, na.rm=TRUE), .groups = "drop") |>
    mutate(n_sub = glue::glue("N: {n_sub}")) |>
    rename(`N Sub` = n_sub) |>
    ggplot(aes(y=contrast, x=d, color=Network)) +
    facet_grid(Task~`N Sub`, shrink = TRUE, scales = "free_y") +
    geom_boxplot(outlier.alpha = 0.25, outlier.size = .5) +
    # geom_jitter(width = 0, alpha=0.1) +
    scale_x_continuous(
      name = "avg dist(Gold Peak, Study Peak) (mm)",
      breaks = c(0, 50),
      labels = c(0, 50),
      limits = c(0, 150)) +
    scale_color_viridis_d(option = "turbo") +
    guides(colour = guide_legend(title.position = "top", nrow=2)) +
    theme_gray(base_size = 10) +
    theme(legend.position = "bottom")
  ggsave("analyses/avg-dist-network-all-p01.png", p, device = ragg::agg_png, width=8, height = 7)

  p2 <- space |>
    filter(fct_match(FWE, "0.01")) |>
    filter(contrast %in% c("MATH-STORY","CUE-AVG","MATCH-REL","RANDOM-TOM","2BK-0BK")) |>
    group_by(Task, contrast, n_sub, iter, Network) |>
    summarise(d = mean(d, na.rm=TRUE), .groups = "drop") |>
    mutate(n_sub = glue::glue("N: {n_sub}")) |>
    rename(`N Sub` = n_sub) |>
    ggplot(aes(y=contrast, x=d, color=Network)) +
    facet_grid(Task~`N Sub`, shrink = TRUE, scales = "free_y") +
    geom_boxplot(outlier.alpha = 0.25, outlier.size = .5) +
    # geom_jitter(width = 0, alpha=0.1) +
    scale_x_continuous(
      name = "avg dist(Gold Peak, Study Peak) (mm)",
      breaks = c(0, 50),
      labels = c(0, 50),
      limits = c(0, 150)) +
    scale_color_viridis_d(option = "turbo") +
    guides(colour = guide_legend(title.position = "top", nrow=3)) +
    theme_gray(base_size = 8) +
    theme(legend.position = "bottom")
  
  ggsave("analyses/avg-dist-network-small.png", p2, device = ragg::agg_png, width=4, height = 6)

  avg <- space |>
    filter(fct_match(FWE, "0.95")) |>
    distinct(x, y, z, Network, Task, contrast, Value) |>
    group_by(Network) |>
    summarise(avg = mean(Value), .groups = "drop")
  
  p3 <- space |>
    filter(fct_match(FWE, "0.95")) |>
    group_by(n_sub, iter, Network) |>
    summarise(
      d = mean(d, na.rm=TRUE), 
      .groups = "drop") |>
    left_join(avg) |>
    rename(`N Sub` = n_sub) |>
    ggplot(aes(y=Network, x=d, color=avg)) +
    facet_wrap(~`N Sub`, ncol=1, labeller = label_both) +
    geom_boxplot(outlier.alpha = 0.25, outlier.size = .5) +
    # geom_jitter(width = 0, alpha=0.1) +
    scale_x_continuous(
      name = "avg dist(Gold Peak, Study Peak) (mm)") +
    scale_color_viridis_c(option = "turbo") +
    theme_gray(base_size = 8) +
    labs(color = "Cohen's d") +
    guides(color=guide_colourbar(title.position = "top", barheight = .5)) +
    theme(legend.position = "bottom")
  ggsave("analyses/avg-dist-network.png", p3, device = ragg::agg_png, width=4, height = 5)
  
}


avg_network_size <- function(){
  store<- "_hcp"
  
  load_space <- function(){
    targets::tar_read(space, store = store) |>
      filter(!is.na(label)) |>
      mutate(corrp_thresh = factor(corrp_thresh)) |>
      rename(Network = `Network Name`) |>
      mutate(
        Value = Value / sqrt(n_pop),
        label = if_else(
          is.na(`Network`), 
          label, 
          str_extract(label, "(?<=_)([[:alnum:]]+)_[[:digit:]]+$")),
        `Network` = if_else(
          is.na(`Network`), 
          "subcortical", 
          `Network`
        ),
        `Network` = factor(
          `Network`, 
          levels = c(
            "limbic",
            "default",
            "control",
            "salience / ventral attention",
            "dorsal attention",
            "somatomotor", "visual",
            "subcortical"
          ),
          ordered = TRUE)
      ) 
  }
  space <- load_space() |>
    left_join(targets::tar_read(contrasts, store=store),by = c("Task", "CopeNumber")) |>
    rename(contrast=ContrastName, FWE = corrp_thresh) 
  
  avg <- space |>
    filter(fct_match(FWE, "0.01")) |>
    distinct(x, y, z, Network, Task, contrast, Value) |>
    group_by(Network, Task, contrast) |>
    summarise(avg = mean(Value), .groups = "drop")
  
  space |>
    filter(fct_match(FWE, "0.01")) |>
    whoppeR::WISEsummary(
      dependentvars = "d",
      betweenvars = c("Network", "n_sub", "Task", "contrast"),
      na.rm = TRUE
    ) |>
    left_join(avg) |>
    rename(`N Sub` = n_sub) |>
    ggplot(aes(y=Network)) +
    facet_grid(Task+contrast~`N Sub`, shrink = TRUE, scales = "free_y") +
    geom_errorbarh(aes(color=avg, xmin=d_CI_lower, xmax=d_CI_upper)) +
    scale_color_viridis_c(option="turbo") +
    scale_x_continuous(
      name = "avg dist(Gold Peak, Study Peak) (mm)") +
    theme_gray(base_size = 10)
}
