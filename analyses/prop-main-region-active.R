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

Sys.setenv(TAR_PROJECT = "hcp")

at <- targets::tar_read(at) |>
  mutate(across(c(x,y,z), as.integer))

at_list <- targets::tar_read(at_list) |>
  mutate(across(c(x,y,z), as.integer))

targets::tar_load(pop_d)

out <- targets::tar_read(prop0all) |>
  dplyr::left_join(at_list) |>
  distinct(
    n_sub, iter, Task, ContrastName, CopeNumber, label, `Label Name`,
    `Network Name`, `Full component name`, n_parcels) |>
  collect() |>
  count(
    n_sub, Task, ContrastName, CopeNumber, label, `Label Name`,
    `Network Name`, `Full component name`, n_parcels) |>
  mutate(prop = n / 100)

# these are the regions whose activation is the "strongest"
# so, assumption is that researcher picks the task to activate
# this specific region
gold <- pop_d |>
  left_join(at_list) |>
  group_by(
    Task, label, `Label Name`,
    `Network Name`, `Full component name`, n_parcels) |>
  summarise(
    d = mean(cope/sigma),
    .groups = "drop"
  ) |>
  group_by(Task, n_parcels) |>
  slice_max(order_by=d, n=10) |>
  ungroup() |>
  filter(!is.na(label)) |> 
  group_by(Task) |>
  mutate(
    l=label |> 
      factor() |> 
      as.numeric() |>
      factor()) |>
  ungroup()

# regions with at least one voxel active
out |>
  right_join(distinct(gold, Task, l, label, n_parcels)) |>
  ggplot(aes(x=n_sub, y=prop, group=l)) +
  geom_point(show.legend = FALSE, alpha=0.2) +
  geom_line(show.legend = FALSE, alpha=0.2) +
  facet_grid(n_parcels~Task) +
  scale_y_continuous(
    "Proportion Simulations w/ Activity in Most Active ROI",
    limits = c(0.25, 1)
  ) +
  xlab("N Sub")

ggsave("prop-active-most-active-roi.png", width = 7, height = 6)


library(ggsegSchaefer)

template <- RNifti::readNifti(MNITemplate::getMNIPath(res="2mm"))
mask <- array(0, dim=dim(template))
tmp <- filter(at,  label=="7Networks_RH_Vis_23")
for (r in 1:nrow(tmp)){
  mask[tmp[r,]$x, tmp[r,]$y, tmp[r,]$z] <- 1
}


quick <- function(task){
  gold  |>
    filter(Task==task) |> 
    rename(region=label) |>
    ggplot() +
    ggseg::geom_brain(
      aes(fill=l),
      position = position_brain(hemi ~ side),
      atlas = ggsegSchaefer::schaefer7_400,
      show.legend = FALSE) +
    scale_fill_viridis_d(option="turbo", na.value = NA) +
    theme_void()
}

a + 
  quick("GAMBLING") + 
  quick("LANGUAGE") + 
  quick("MOTOR") + 
  quick("RELATIONAL") + 
  quick("SOCIAL") + 
  quick("WM") + 
  plot_layout(
    design = c(
      area(t = 1, l = 1, b = 11, r = 6),
      area(t = 3, l = 1, b = 5, r = 2),
      area(t = 3, l = 3, b = 5, r = 4),
      area(t = 3, l = 5, b = 5, r = 6),
      area(t = 9, l = 1, b = 11, r = 2),
      area(t = 9, l = 3, b = 11, r = 4),
      area(t = 9, l = 5, b = 11, r = 6)
    )
  )

ggsave("prop-top-regions.png")


targets::tar_load(space)



# space |>
#   filter(
#     corrp_thresh==0.95,
#     !is.na(study_ind)) |> 
#   semi_join(gold_peaks) |>
#   left_join(
#     at |> 
#       mutate(`Label Name` = if_else(is.na(`Label Name`), label, `Label Name`)) |>
#       select(x, y, z, label.study=`Label Name`),
#     by = c("x.study"="x", "y.study"="y", "z.study"="z")) |>
#   mutate(
#     `Label Name` = if_else(is.na(`Label Name`), "White Matter", `Label Name`),
#     label.study = if_else(is.na(label.study), "White Matter",label.study )) |>
#   group_by(Task, n_sub, x, y, z, iter, `Label Name`) |>
#   summarize(
#     any_simulations = any(`Label Name` == label.study),
#     .groups = "drop"
#   ) |>
#   group_by(Task, n_sub, x, y, z) |>
#   summarize(any_simulations=sum(any_simulations), .groups="drop") |>
#   unite(col="peak", x, y, z) |>
#   mutate(any_simulations = any_simulations / 100) |>
#   ggplot(aes(x=n_sub, group=peak, y=any_simulations)) +
#   geom_point(alpha=0.2) +
#   geom_line(alpha=0.2) +
#   facet_wrap(~Task) +
#   scale_y_continuous(
#     "% Studies with Peak",
#     limits = c(0,1),
#     breaks = c(0,0.5,1),
#     labels = c(0,0.5,1)
#   )


space |>
  filter(corrp_thresh==0.95) |>
  select(Task, iter, n_sub, x=x.study, y=y.study, z=z.study) |>
  left_join(select(at_list, -`Label Name`, -hemi, -n_voxels, -volume)) |>
  semi_join(select(gold, Task, label, n_parcels)) |>
  distinct(Task, iter, label, n_sub, n_parcels) |>
  right_join(select(gold, Task, label, n_parcels) |> crossing(distinct(space, n_sub))) |>
  group_by(Task, label, n_sub, n_parcels) |>
  summarise(
    prop = sum(!is.na(iter)) / 100,
    .groups = "drop"
  ) |>
  ggplot(aes(x=n_sub, y=prop, group=label)) +
  geom_point(show.legend = FALSE, alpha=0.2) +
  geom_line(show.legend = FALSE, alpha=0.2) +
  facet_grid(n_parcels~Task) +
  scale_y_continuous(
    "Proportion Simulations w/ Peak in Most Active ROI",
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c(0, 1)
  ) +
  xlab("N Sub")

ggsave("prop-w-peak-most-active-roi", width = 7, height = 6)


gold2 <- targets::tar_read(gold_peaks) |>
  select(Task, m) |>
  unnest(m) |>
  left_join(at_list) |>
  filter(!is.na(label)) |> 
  group_by(Task, label, n_parcels, n_networks) |>
  slice_max(order_by = Value, n = 1) |>
  group_by(Task, n_parcels, n_networks) |>
  slice_max(order_by = Value, n = 10, with_ties = FALSE) |>
  ungroup() |>
  group_by(Task) |>
  mutate(
    l=label |> 
      factor() |> 
      as.numeric() |>
      factor()) |>
  ungroup()

  
space |>
  filter(corrp_thresh==0.95) |>
  select(Task, iter, n_sub, x=x.study, y=y.study, z=z.study) |>
  left_join(select(at_list, -`Label Name`, -hemi, -n_voxels, -volume)) |>
  semi_join(distinct(gold2, Task, label, n_parcels)) |>
  distinct(Task, iter, label, n_sub, n_parcels) |>
  right_join(distinct(gold2, Task, label, n_parcels) |> crossing(distinct(space, n_sub))) |>
  group_by(Task, label, n_sub, n_parcels) |>
  summarise(
    prop = sum(!is.na(iter)) / 100,
    .groups = "drop"
  ) |>
  ggplot(aes(x=n_sub, y=prop, group=label)) +
  geom_point(show.legend = FALSE, alpha=0.2) +
  geom_line(show.legend = FALSE, alpha=0.2) +
  facet_grid(n_parcels~Task) +
  scale_y_continuous(
    "Proportion Simulations with Peak in ROI with Peak",
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c(0, 1)
  ) +
  xlab("N Sub")

ggsave("prop-w-peak-roi-w-peak.png", width = 7, height = 6)

# for the most active regions, which have the highest peaks
gold_peaks0 <- targets::tar_read(gold_peaks) |>
  select(Task, m) |>
  unnest(m) |>
  left_join(at_list) |>
  filter(!is.na(label)) |> 
  group_by(Task, label, n_parcels, n_networks) |>
  slice_max(order_by = Value, n = 1, with_ties = FALSE) |>
  ungroup() |>
  distinct(Task, n_parcels, n_networks, label) |>
  mutate(local_peak = TRUE)


golden2 <- pop_d |>
  left_join(at_list) |>
  group_by(
    Task, label, `Label Name`,
    `Network Name`, `Full component name`, n_parcels, n_networks) |>
  summarise(
    d = mean(cope/sigma),
    peak = max(cope / sigma),
    .groups = "drop"
  ) |>
  filter(!is.na(label)) |>
  left_join(gold_peaks0) |>
  replace_na(list(local_peak = FALSE)) |>
  rename(`Local Peak` = local_peak)


# double check that there are the correct number of regions
golden2 |>
  count(Task, n_parcels)

golden2 |>
  filter(d > 0) |>
  count(Task, n_parcels)

maxes <- golden2 |>
  filter(n_parcels==400) |>
  group_by(Task, n_parcels) |>
  slice_max(d, n=10, with_ties = FALSE) |>
  slice_min(d, n=1, with_ties = FALSE) |>
  ungroup() |>
  distinct(Task, d, n_parcels)|>
  rename(top=d)

a <- golden2 |>
  filter(d>0) |>
  ggplot(aes(x=d, y=peak, color=`Local Peak`)) +
  geom_point(alpha=0.1) +
  facet_grid( n_parcels ~ Task) +
  xlab("Average d") +
  ylab("Maximum d") +
  scale_x_continuous(
    n.breaks = 3,
    limits = c(0, 2)
  ) +
  scale_y_continuous(
    limits = c(0, 2),
    n.breaks = 3
  )

b <- golden2 |>
  filter(d>0) |>
  pivot_longer(c(d, peak)) |>
  mutate(
    name = case_when(
      name == "d" ~ "average",
      name == "peak" ~ "max"
    )) |>
  ggplot(aes(color=`Local Peak`, x=name, y=value)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
  facet_grid( n_parcels ~ Task) +
  xlab("Statistic") +
  ylab("Cohen's d") +
  scale_y_continuous(
    limits = c(0, 2),
    n.breaks = 3
  )

a + b + plot_layout(guides = "collect") & theme(legend.position="bottom")
ggsave("avg-max.png", width=12, height = 8)


# question: how can a region have a maximum without having a local peak?
# answer: voxels are always on the border!
select(golden2, -`Label Name`, -`Network Name`) |> 
  filter(!`Local Peak`, n_parcels==400, Task=="MOTOR")

tstat_tbl <- to_tbl(
  "/home/ubuntu/mnt/meta/meta/data-raw/hcp-niis/nsub-480_iter-0_flags-AVG_tstat1.nii.gz"
  ) |> 
  left_join(at_list) |> 
  filter(!is.na(label), n_parcels==400) |>
  mutate(x=x-1, y=y-1, z=z-1)

tstat_tbl |>
  filter(label=="7Networks_LH_Cont_Cing_1")

pop_d |> 
  mutate(x=x-1, y=y-1, z=z-1) |>
  filter(Task=="MOTOR", x==47, y==73, z==46)

tstat_tbl |>
  group_by(label) |>
  slice_max(order_by = value, n=1) |>
  ungroup() |>
  filter(label=="7Networks_LH_Cont_Par_1")
