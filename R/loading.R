

# load_tab <- function(fname, fields, n_max){
#   readr::read_tsv(
#     fname,
# col_select = c(
#   "f.eid",
#   tidyselect::starts_with(stringr::str_c("f.", fields))),
#     n_max = n_max) |>
#     dplyr::distinct(f.eid, .keep_all = TRUE) |>
#     tidyr::pivot_longer(
#       cols = c(tidyselect::starts_with("f."), -"f.eid"),
#       names_to = c("field", "instance", "array"),
#       names_pattern = "f.([[:digit:]]+).([[:digit:]]).([[:digit:]])",
#       names_transform = list(field = as.double, instance=as.integer, array=as.integer)) |>
#     tidyr::pivot_wider(names_from = field) |>
#     na.omit()
#   
# }

load_key <- function(fname){
  readr::read_csv(
    fname,
    col_select = c("FieldID", "Field", "ValueType", "Units", "Category"),
    col_types = readr::cols(
      FieldID = readr::col_integer())) |>
    dplyr::filter(
      (Category %in% c(110, 112, 119, 106, 107, 111, 109)) | (FieldID %in% c(20016, 21022, 21000)),
      !stringr::str_detect(Field, "NIFTI"),
      !stringr::str_detect(Field, "eprime|Eprime"),
      !stringr::str_detect(Field, "position"),
      !stringr::str_detect(Field, "correlation"),
      !stringr::str_detect(Field, "BOLD effect"),
      !stringr::str_detect(Field, "head motion"),
      !stringr::str_detect(Field, "percentile"),
      !stringr::str_detect(Field, "Discrepancy"),
      !stringr::str_detect(Field, "Amount of warping"),
      !stringr::str_detect(Field, "scaling from T1 head"),
      !stringr::str_detect(Field, "T2star"),
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
}

load_tab <- function(fname, key, n_max){
  readr::read_tsv(
    fname,
    n_max = n_max,
    col_select = c(
      "f.eid",
      tidyselect::matches(glue::glue("f.({paste(unique(key$FieldID), sep='|')}).[[:digit:]]+.[[:digit:]]+"))),  
    col_types = readr::cols(
      f.eid = readr::col_integer(),
      .default = readr::col_character())) |>
    dplyr::filter(!(is.na(f.20016.0.0) & is.na(f.20016.1.0) & is.na(f.20016.2.0) & is.na(f.20016.3.0))) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      f.20016.0.0 = max(c(f.20016.0.0, f.20016.1.0, f.20016.2.0, f.20016.3.0), na.rm=TRUE)) |>
    dplyr::ungroup() |>
    dplyr::select(-f.21000.1.0, -f.21000.2.0, -f.20016.1.0, -f.20016.2.0, -f.20016.3.0) |>
    na.omit() |>
    dplyr::rename(age = `f.21022.0.0`, ethnicity = `f.21000.0.0`) |>
    tidyr::pivot_longer(
      cols = c(tidyselect::starts_with("f."), -"f.eid", -age, -ethnicity),
      names_to = c("FieldID", "instance", "array"),
      names_pattern = "f.([[:digit:]]+).([[:digit:]]+).([[:digit:]]+)",
      names_transform = list(FieldID = as.integer, instance=as.integer, array=as.integer)) |>
    dplyr::select(-instance, -array) |>
    dplyr::semi_join(key, by = "FieldID") |>
    tidyr::pivot_wider(names_from = FieldID) |>
    dplyr::mutate(
      dplyr::across(.cols = tidyselect::matches("[[:digit:]]+"), .fns = as.numeric),
      age = as.numeric(age),
      ethnicity = as.numeric(ethnicity)) |>
    dplyr::mutate(
      ethnicity = dplyr::case_when(
        ethnicity %in% c(1, 1001, 1002, 1003) ~ "white",
        ethnicity %in% c(2, 2003, 2004) ~ "multiple",
        ethnicity %in% c(4001) ~ "caribbean",
        ethnicity %in% c(2001, 2002, 4002, 4003, 4) ~ "black",
        ethnicity %in% c(3001) ~ "indian",
        ethnicity %in% c(5) ~ "chinese",
        ethnicity %in% c(3, 3004, 3003, 3002) ~ "other asian",
        TRUE ~ "other"
      ),
      age = cut(age, seq(30, 70, by=5))) |>
    dplyr::distinct(f.eid, ethnicity, .keep_all = TRUE)
  
}