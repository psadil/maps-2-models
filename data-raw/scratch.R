library(tidyverse)

key <- readr::read_csv(
  "data-raw/Data_Dictionary_Showcase.csv",
  col_select = c("FieldID", "Field", "ValueType", "Units", "Category"),
  col_types = cols(
    FieldID = col_integer())) |>
  filter(
    (Category %in% c(110, 112, 119, 106, 107, 111, 109)) | (FieldID == 20016),
    !str_detect(Field, "NIFTI"),
    !str_detect(Field, "eprime|Eprime"),
    !str_detect(Field, "position"),
    !str_detect(Field, "correlation"),
    !str_detect(Field, "BOLD effect"),
    !str_detect(Field, "head motion"),
    !str_detect(Field, "percentile"),
    !str_detect(Field, "Discrepancy"),
    !str_detect(Field, "Amount of warping"),
    !str_detect(Field, "scaling from T1 head"),
    !str_detect(Field, "T2star"),
    !str_detect(Field, "activation"),
    !str_detect(Field, "omponent amplitudes"))

# for list of imaging categories: https://biobank.ctsu.ox.ac.uk/crystal/label.cgi?id=100

d <- readr::read_csv(
  here::here("data-raw/ukb37917_narrow.csv"), 
  n_max = 1000,
  # col_select = c(
  #   "f.eid",
  #   tidyselect::matches(glue::glue("f.({paste(unique(key$FieldID), sep='|')}).[[:digit:]]+.[[:digit:]]+"))),  
  col_types = cols(
    f.eid = col_integer(),
    .default = col_character())) |>
  dplyr::distinct(f.eid, .keep_all = TRUE) |>
  tidyr::pivot_longer(
    cols = c(tidyselect::starts_with("f."), -"f.eid"),
    names_to = c("FieldID", "instance", "array"),
    names_pattern = "f.([[:digit:]]+).([[:digit:]]+).([[:digit:]]+)",
    names_transform = list(FieldID = as.integer, instance=as.integer, array=as.integer)) |>
  na.omit() |>
  semi_join(key, by = "FieldID") |>
  tidyr::pivot_wider(names_from = FieldID) |>
  na.omit() |>
  mutate(across(.cols = matches("[[:digit:]]+"), .fns = as.numeric))


