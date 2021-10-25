library(tidyverse)

# source(here::here("R", "loading.R"))

key <- readr::read_csv(
  here::here("data-raw","Data_Dictionary_Showcase.csv"),
  col_select = c("FieldID", "Field", "ValueType", "Units", "Category"),
  col_types = readr::cols(
    FieldID = readr::col_integer())) %>%
  dplyr::filter(
    (Category %in% c(119, 106, 111, 109)) | 
      (FieldID %in% c(20016, 21022, 21000, 20126, 25010, 20252, 23104, 25741, 21001, 130874, 42020)),
    !(stringr::str_detect(Field, "NIFTI") & !str_detect(Field, "T1")),
    !stringr::str_detect(Field, "eprime|Eprime"),
    !stringr::str_detect(Field, "position"),
    !stringr::str_detect(Field, "correlation"),
    !stringr::str_detect(Field, "BOLD effect"),
    !stringr::str_detect(Field, "percentile"),
    !stringr::str_detect(Field, "Discrepancy"),
    !stringr::str_detect(Field, "Amount of warping"),
    !stringr::str_detect(Field, "scaling from T1 head"),
    !stringr::str_detect(Field, "T2star"),
    !stringr::str_detect(Field, "T2-FLAIR"),
    !stringr::str_detect(Field, "activation"),
    !stringr::str_detect(Field, "omponent amplitudes"),
    !stringr::str_detect(Field, "Intensity scaling"),
    !stringr::str_detect(Field, "Increased search space"),
    !stringr::str_detect(Field, "outlier slices"),
    !stringr::str_detect(Field, "-to-noise ratio"),
    !stringr::str_detect(Field, "Echo Time"),
    !stringr::str_detect(Field, "faces-shapes"),
    !stringr::str_detect(Field, "DICOM"),
    !stringr::str_detect(Field, "surface model files"))

# "/dcl01/smart/data/UKBiobank/BehavioralData/ukb41898.tab",
# 
# "data-raw/ukb37917.tab",
# "data-raw/ukb26883.tab",
# n_max = 1000,

d <- readr::read_tsv(
  num_threads = 1,
  # "data-raw/ukb47934.tab",
  # "data-raw/ukb26883.tab",
  # "data-raw/ukb37917.tab",
  # n_max = 300000,
  # "/dcl01/smart/data/UKBiobank/BehavioralData/ukb47934.tab",
  "/dcl01/smart/data/UKBiobank/BehavioralData/ukb37917.tab",
  # "/dcl01/smart/data/UKBiobank/BehavioralData/ukb26883.tab",
  col_select = c(
    "f.eid",
    tidyselect::matches(glue::glue("f.({paste(unique(key$FieldID), sep='|')}).[[:digit:]]+.[[:digit:]]+"))),  
  col_types = readr::cols(
    f.eid = readr::col_integer(),
    .default = readr::col_character())) %>%
  tidyr::pivot_longer(
    cols = c(tidyselect::starts_with("f."), -"f.eid"),
    names_to = c("FieldID", "instance", "array"),
    names_pattern = "f.([[:digit:]]+).([[:digit:]]+).([[:digit:]]+)",
    names_transform = list(FieldID = as.integer, instance=as.integer, array=as.integer)) %>%
  na.omit() %>%
  group_by(f.eid, FieldID) %>%
  summarise(value = max(value), .groups = "drop") %>%
  dplyr::semi_join(key, by = "FieldID") %>%
  tidyr::pivot_wider(names_from = FieldID) %>%
  rename(age = `21022`, ethnicity = `21000`) 
# distinct(f.eid, ethnicity, .keep_all = TRUE)


d <- read_csv("data-raw/ukb37917_narrow.csv") %>%
  mutate(
    `20252` = !is.na(`20252`),
    ethnicity = dplyr::case_when(
      ethnicity %in% c(1, 1001, 1002, 1003) ~ "white",
      ethnicity %in% c(2, 2003, 2004) ~ "multiple",
      ethnicity %in% c(4001) ~ "caribbean",
      ethnicity %in% c(2001, 2002, 4002, 4003, 4) ~ "black",
      ethnicity %in% c(3001) ~ "indian",
      ethnicity %in% c(5) ~ "chinese",
      ethnicity %in% c(3, 3004, 3003, 3002) ~ "other asian",
      TRUE ~ "other"),
    `20126` = factor(
      `20126`,
      labels = c(
        "nothing",
        "Bipolar 1",
        "Bipolar 2",
        "recurrent depression (severe)",
        "recurrent depression (moderate)",
        "isolated depression")),
    age = cut(age, seq(30, 70, by=5))) %>%
  rename(T1w = `20252`, depression = `20126`)


d |> 
  ggplot(aes(x=ethnicity)) + 
  facet_wrap(~T1w, labeller = label_both, scales = "free") +
  geom_bar()

d %>%
  filter(T1w) %>%
  group_by(ethnicity) %>%
  summarise(N = n())

d %>%
  filter(T1w) %>%
  group_by(depression) %>%
  summarise(N = n(), .groups = "drop")


d |> 
  ggplot(aes(x=depression)) + 
  geom_bar() + 
  facet_wrap(~ethnicity, scales = "free_x") +
  coord_flip()

d |> 
  filter(T1w) %>%
  ggplot(aes(x=depression)) + 
  geom_bar() + 
  facet_wrap(~ethnicity, scales = "free_x") +
  coord_flip()


d %>% 
  filter(T1w, !is.na(depression)) %>%
  ggplot() +
  facet_wrap(~depression, scales = "free") +
  geom_bar(aes(x=age)) 


d |>
  ggplot() +
  geom_bar(aes(x=age)) +
  facet_wrap(~ ethnicity, scales = "free") 


tar_load(ukb)

# brain size, normalized for head size
ukb |> 
  pivot_longer(c(-f.eid, -ethnicity, -age, -`20016`)) |> 
  filter(name=="25009") |> 
  ggplot(aes(x=value, y=`20016`)) + 
  geom_point(alpha = 0.05) + 
  facet_wrap(~ethnicity) + 
  stat_smooth(method = "lm")

# brain size, not normalized
ukb |> 
  pivot_longer(c(-f.eid, -ethnicity, -age, -`20016`)) |> 
  filter(name=="25010") |> 
  ggplot(aes(x=value, y=`20016`)) + 
  geom_point(alpha = 0.05) + 
  facet_wrap(~ethnicity) + 
  stat_smooth(method = "lm") +
  xlab("Brain Volume (mm^3)") +
  scale_y_continuous(
    name = "Fluid Intelligence (0-14)",
    limits = c(0, 14)) 


p <- ukb |> 
  # slice_sample(n = 1000) |>
  pivot_longer(c(-f.eid, -ethnicity, -age, -`20016`)) |> 
  filter(name=="25010") |> 
  # mutate(ethnicity = fct_drop(ethnicity)) |>
  ggstatsplot::grouped_ggscatterstats(
    x = value,
    y = `20016`,
    type = "Non-parametric",
    grouping.var = ethnicity, # grouping variable
    # label.expression = length > 200,
    xlab = "brain volume (mm^3)",
    ylab = "fluid intelligence (0-14)",
    ggtheme = ggplot2::theme_grey(base_size = 8),
    ggplot.component = list(
      ggplot2::scale_y_continuous(
        limits = c(0, 14)),
      ggplot2::scale_x_continuous(
        limits = c(820000, 1700000)))
    
    # plotgrid.args = list(nrow = 1),
    # annotation.args = list(title = "Relationship between movie length and IMDB ratings")
  ) 

ggsave("fi_by-ethnicity.png", p, device = ragg::agg_png(), width = 20, height = 9, units = "in")



# depression

d <- readr::read_tsv(
  "data-raw/ukb37917.tab",
  n_max = 50000,
  col_select = c(
    "f.eid",
    "f.21022.0.0", 
    "f.21000.0.0",
    "f.20126.0.0"),  
  col_types = readr::cols(
    .default = readr::col_integer())) |>
  tidyr::pivot_longer(
    cols = c(tidyselect::starts_with("f."), -"f.eid"),
    names_to = c("FieldID", "instance", "array"),
    names_pattern = "f.([[:digit:]]+).([[:digit:]]+).([[:digit:]]+)",
    names_transform = list(FieldID = as.integer, instance=as.integer, array=as.integer)) |>
  na.omit() |>  
  pivot_wider(names_from = "FieldID") |>
  rename(age = `21022`, ethnicity = `21000`, depression = `20126`) |>
  mutate(
    ethnicity = dplyr::case_when(
      ethnicity %in% c(1, 1001, 1002, 1003) ~ "white",
      ethnicity %in% c(2, 2003, 2004) ~ "multiple",
      ethnicity %in% c(4001) ~ "caribbean",
      ethnicity %in% c(2001, 2002, 4002, 4003, 4) ~ "black",
      ethnicity %in% c(3001) ~ "indian",
      ethnicity %in% c(5) ~ "chinese",
      ethnicity %in% c(3, 3004, 3003, 3002) ~ "other asian",
      TRUE ~ "other"),
    age = cut(age, seq(30, 70, by=5)),
    depression = 
      factor(
        depression,
        labels = c(
          "nothing",
          "Bipolar 1",
          "Bipolar 2",
          "recurrent depression (severe)",
          "recurrent depression (moderate)",
          "isolated depression")))

