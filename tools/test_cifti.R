

roi_from_cifti <- function(file, name="Schaefer_400"){
  parc <- ciftiTools::load_parc(name)
  parc <- ciftiTools::parc_add_subcortex(parc)
  xii <- ciftiTools::read_xifti(file)
  ciftiTools::apply_parc(xii, parc, FUN=mean, na.rm=TRUE) |>
    tibble::as_tibble(rownames = "label") |>
    dplyr::rename(Z = V1) |>
    dplyr::filter(!label == "???") |>
    na.omit()
}


tbl_from_cifti <- function(file){
  ciftiTools::read_xifti(file, flat = TRUE) |>
    t() |>
    dplyr::as_tibble() |>
    dplyr::mutate(t = dplyr::row_number()) |>
    tidyr::pivot_longer(-t, names_prefix = "V") |>
    dplyr::mutate(name = as.integer(name))
}




