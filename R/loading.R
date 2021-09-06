

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
      (Category %in% c(110, 112, 119, 106, 107, 111, 109)) | (FieldID == 20016),
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
    dplyr::distinct(f.eid, .keep_all = TRUE) |>
    tidyr::pivot_longer(
      cols = c(tidyselect::starts_with("f."), -"f.eid"),
      names_to = c("FieldID", "instance", "array"),
      names_pattern = "f.([[:digit:]]+).([[:digit:]]+).([[:digit:]]+)",
      names_transform = list(FieldID = as.integer, instance=as.integer, array=as.integer)) |>
    dplyr::semi_join(key, by = "FieldID") |>
    tidyr::pivot_wider(names_from = FieldID) |>
    na.omit() |>
    dplyr::mutate(dplyr::across(.cols = matches("[[:digit:]]+"), .fns = as.numeric))
  
}