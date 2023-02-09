Sys.setenv(TAR_PROJECT = "tfce")
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
source(here::here("R","utils.R"))
source(here::here("R","tfce.R"))
source(here::here("R","updates.R"))
source(here::here("R","poster.R"))

n_pop <- targets::tar_read(tfce_pop)$n_sub

# mri_vol2vol --targ $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz --mov ~/mnt/Yeo_JNeurophysiol11_MNI152/Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz --o ~/mnt/yeo_MNI152_1mm.nii.gz --regheader --no-save-reg

# yeo <- to_tbl("/home/ubuntu/mnt/yeo7_MNI152_2mm.nii.gz", measure = "Network Order") |>
#   left_join(readr::read_csv("data-raw/1000subjects_reference/7NetworksOrderedNames.csv"))
# yeo_labels <- readr::read_csv("data-raw/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/Yeo2011_7networks_N1000.split_components.glossary.csv") |>
#   mutate(value = 1:n()) |>
#   right_join(to_tbl("data-raw/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/MNI152/Yeo2011_7Networks_N1000.split_components.FSL_MNI152_2mm.nii.gz") )

space <- targets::tar_read(space) |>
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
        "limbic A",
        "default A", "default B", "default C",
        "control A", "control B",
        "salience A", "salience B",
        "dorsal atten A", "dorsal atten B",
        "temporal parietal",
        "somatomotor A", "somatomotor B", "central vis", "peripheral vis",
        "subcortical"
      ),
      ordered = TRUE)) 


# space |>
#   filter(!is.na(`Network Name`)) |>
#   ggplot(aes(x=nsub, y=d)) +
#   facet_wrap(~`Network Name`, labeller = labeller(nm = label_wrap_gen(width = 25))) +
#   geom_hline(yintercept =0, color="darkgoldenrod1", size=1) +
#   geom_boxplot(outlier.alpha = 0.5) +
#   ylab("dist(Gold Peak, Study Peak) (mm)") +
#   xlab("N Participants Simulated")


space  |>
  ggplot(aes(y=`Network`, x=d)) +
  facet_wrap(~n_sub,nrow = 1) +
  geom_boxplot(outlier.alpha = 0.5)

space |>
  ggplot(aes(y=label, x=d)) +
  facet_grid(
    `Network`~n_sub, 
    scales = "free_y",space = "free_y", 
    labeller = labeller(`Network` = label_wrap_gen(width = 13))) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme_gray(base_size = 8)


space |>
  whoppeR::WISEsummary(
    dependentvars = "d", 
    betweenvars = c("hemi", "Network", "n_sub")) |>
  ggplot(aes(x=d_mean, y=Network, color=hemi)) +
  facet_wrap(~n_sub,nrow=1) +
  geom_errorbarh(aes(xmin=d_CI_lower, xmax=d_CI_upper))


space |>
  whoppeR::WISEsummary(
    dependentvars = "d", 
    betweenvars = c("hemi", "label", "Network", "n_sub")) |>
  ggplot(aes(x=d_mean, y=label, color=hemi)) +
  facet_grid(
    `Network`~n_sub, 
    scales = "free_y",space = "free_y", 
    labeller = labeller(`Network` = label_wrap_gen(width = 13))) +
  geom_errorbarh(aes(xmin=d_CI_lower, xmax=d_CI_upper))



space |>
  whoppeR::WISEsummary(
    dependentvars = "d", 
    betweenvars = c("label", "Value", "Network", "n_sub", "hemi")) |>
  ggplot(aes(x=d_mean, y=label, color=Value, shape=hemi)) +
  geom_point() +
  facet_grid(
    `Network`~n_sub, 
    scales = "free_y",space = "free_y", 
    labeller = labeller(`Network` = label_wrap_gen(width = 13))) +
  geom_errorbarh(aes(xmin=d_CI_lower, xmax=d_CI_upper))  +
  scale_color_viridis_c(option="turbo")



