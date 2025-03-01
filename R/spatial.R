make_atlas <- function(labelxml, nii) {
  labels <- xml2::read_xml(labelxml) |>
    xml2::as_list() |>
    magrittr::extract2("atlas") |>
    magrittr::extract2("data") |>
    tibble::as_tibble(.name_repair = "unique") |>
    tidyr::pivot_longer(cols = tidyselect::everything()) |>
    dplyr::mutate(
      index = stringr::str_extract(name, "[[:digit:]]+") |> as.integer(),
      label = purrr::map_chr(value, ~.x)
    ) |>
    dplyr::select(index, label)
  
  neurobase::readnii(nii) |>
    neurobase::img_indices(add_values = TRUE) |>
    tibble::as_tibble() |>
    dplyr::rename(index = value) |>
    dplyr::filter(index > 0) |>
    dplyr::left_join(labels, by = "index")
}


rename_labels <- function(.data){
  .data |>
    dplyr::mutate(
      label = dplyr::case_match(
        label,
        "Right Accumbens" ~ "Accumbens-R",
        "Left Accumbens" ~ "Accumbens-L",
        "Right Amygdala" ~ "Amygdala-R",
        "Left Amygdala" ~ "Amygdala-L",
        "Right Putamen" ~ "Putamen-R",
        "Left Putamen" ~ "Putamen-L",
        "Right Pallidum" ~ "Pallidum-R",
        "Left Pallidum" ~ "Pallidum-L",
        "Right Hippocampus" ~ "Hippocampus-R",
        "Left Hippocampus" ~ "Hippocampus-L",
        "Right Caudate" ~ "Caudate-R",
        "Left Caudate" ~ "Caudate-L",
        "Right Thalamus" ~ "Thalamus-R",
        "Left Thalamus" ~ "Thalamus-L",
        "Brain-Stem" ~ "Brain Stem",
        .default = label
      )
    ) |>
    dplyr::filter(stringr::str_detect(label, "Cerebell", negate=TRUE)) |>
    dplyr::filter(stringr::str_detect(label, "Diencephalon", negate=TRUE))
}


make_atlas_full <- function(networks = 7, n_parcels = 400) {
  # cortical labels for Schaefer2018_400Parcels_17Networks_order.lut
  # subcortical from harvard-oxford
  
  # always remove all brain-setm/white matter/cortex from subcortical, and always
  # defer to cortical when available
  c_labels <- readr::read_delim(
    here::here(
      "data-raw", "Parcellations", "MNI", "fsleyes_lut",
      glue::glue("Schaefer2018_{n_parcels}Parcels_{networks}Networks_order.lut")
    ),
    col_names = c("index", "R", "G", "B", "label"),
    show_col_types = FALSE
  ) |>
    dplyr::select(index, label) |>
    dplyr::mutate(
      `Label Name` = stringr::str_remove(label, "_[[:digit:]]+$")
    ) |>
    dplyr::left_join(
      readr::read_csv(
        here::here(
          "data-raw",
          "1000subjects_reference",
          "Yeo_JNeurophysiol11_SplitLabels",
          glue::glue("Yeo2011_{networks}networks_N1000.split_components.glossary.csv")
        ),
        show_col_types = FALSE
      ),
      by = "Label Name"
    )
  
  cortical <- to_tbl(
    here::here(
      "data-raw",
      "Parcellations",
      "MNI",
      glue::glue("Schaefer2018_{n_parcels}Parcels_{networks}Networks_order_FSLMNI152_2mm.nii.gz")
    ),
    "index"
  ) |>
    dplyr::filter(index > 0) |>
    dplyr::left_join(c_labels, by = "index")
  
  subcortical <- make_atlas(
    fs::path(fslr::fsldir(), "data/atlases/HarvardOxford-Subcortical.xml"),
    fs::path(fslr::fsldir(), "data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr0-2mm.nii.gz")
  ) |>
    dplyr::filter(stringr::str_detect(label, "Cortex|Matter|Ventricle", negate = TRUE)) |>
    dplyr::anti_join(cortical, by = c("x", "y", "z")) # take only voxels not otherwise accounted for
  
  
  dplyr::bind_rows(cortical, subcortical) |>
    dplyr::mutate(
      hemi = dplyr::if_else(x > 45, "Left", "Right")
    ) |>
    dplyr::select(-index) |>
    dplyr::group_by(label, hemi) |>
    dplyr::mutate(
      n_voxels = dplyr::n(),
      volume = n_voxels * 2^3
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      n_parcels = n_parcels,
      n_networks = networks
    ) |>
    rename_labels()
}

make_atlas_full2 <- function(networks = 7, n_parcels = 400) {
  # cortical labels for Schaefer2018_400Parcels_17Networks_order.lut
  # subcortical from harvard-oxford
  
  # always remove all brain-setm/white matter/cortex from subcortical, and always
  # defer to cortical when available
  c_labels <- readr::read_delim(
    here::here(
      "data-raw", "Parcellations", "MNI", "fsleyes_lut",
      glue::glue("Schaefer2018_{n_parcels}Parcels_{networks}Networks_order.lut")
    ),
    col_names = c("index", "R", "G", "B", "label"),
    show_col_types = FALSE
  ) |>
    dplyr::select(index, label) |>
    dplyr::mutate(
      `Label Name` = stringr::str_remove(label, "_[[:digit:]]+$")
    ) |>
    dplyr::left_join(
      readr::read_csv(
        here::here(
          "data-raw",
          "1000subjects_reference",
          "Yeo_JNeurophysiol11_SplitLabels",
          glue::glue("Yeo2011_{networks}networks_N1000.split_components.glossary.csv")
        ),
        show_col_types = FALSE
      ),
      by = "Label Name"
    )
  
  cortical <- to_tbl(
    here::here(
      "data-raw",
      "Parcellations",
      "MNI",
      glue::glue("Schaefer2018_{n_parcels}Parcels_{networks}Networks_order_FSLMNI152_2mm.nii.gz")
    ),
    "index"
  ) |>
    dplyr::filter(index > 0) |>
    dplyr::left_join(c_labels, by = "index")
  
  subcortical <- make_atlas(
    fs::path(fslr::fsldir(), "data/atlases/HarvardOxford-Subcortical.xml"),
    fs::path(fslr::fsldir(), "data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr0-2mm.nii.gz")
  ) |>
    dplyr::filter(stringr::str_detect(label, "Brain-Stem|Cortex|Matter|Ventricle", negate = TRUE)) |>
    dplyr::mutate(label = stringr::str_remove(label, "Left |Right ")) |>
    dplyr::anti_join(cortical, by = c("x", "y", "z")) # take only voxels not otherwise accounted for
  
  dplyr::bind_rows(cortical, subcortical) |>
    dplyr::mutate(
      hemi = dplyr::if_else(x > 45, "Left", "Right")
    ) |>
    dplyr::group_by(label, hemi) |>
    dplyr::mutate(
      n_voxels = dplyr::n(),
      volume = n_voxels * 2^3
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      n_parcels = n_parcels,
      n_networks = networks
    ) |>
    rename_labels()
}

get_sizes <- function(cls) {
  sizes <- neurobase::img_indices(cls$osize, add_values = TRUE) |>
    tibble::as_tibble() |>
    dplyr::rename(cluster_size = value)
  neurobase::img_indices(cls$oindex, add_values = TRUE) |>
    tibble::as_tibble() |>
    dplyr::filter(value > 0) |>
    dplyr::rename(`Cluster Index` = value) |>
    dplyr::left_join(sizes, by = c("x", "y", "z")) |>
    dplyr::distinct(`Cluster Index`, cluster_size)
}


augment_distance <- function(study, reference, vox_mm = 2) {
  # happens when no peaks were found in study
  if (nrow(study) == 0) {
    return(reference)
  }
  reference |>
    dplyr::rowwise() |>
    dplyr::mutate(study_ind = which.min(sqrt((x - study$x)^2 + (y - study$y)^2 + (z - study$z)^2))) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      x.study = study$x[study_ind],
      y.study = study$y[study_ind],
      z.study = study$z[study_ind],
      d = vox_mm * sqrt((x - x.study)^2 + (y - y.study)^2 + (z - z.study)^2)
    )
}

augment_distance2 <- function(maxes, gold_peaks){
  
  gold <- dplyr::group_nest(gold_peaks, Task)
  
  maxes |>
    dplyr::group_nest(Task, n_sub, iter, fwe_correction) |>
    dplyr::left_join(
      gold, 
      by = dplyr::join_by(Task), 
      suffix = c(".study", ".ref")
    ) |>
    dplyr::mutate(
      m = purrr::map2(
        data.study, data.ref,
        ~augment_distance(
          study=.x,
          reference=.y) |>
          dplyr::select(-n_sub, -iter)
      )
    ) |>
    dplyr::select(-tidyselect::starts_with("data")) |>
    tidyr::unnest(m)
}

add_labels <- function(space, at = make_atlas_full()) {
  dplyr::left_join(space, at, by = c("x", "y", "z"))
  # study peaks can be outside gray matter
}



load_parc <- function(n_parcels=400, src="data-raw/parc.rda") {
  
  # Load the parcellation
  load(src)
  p <- parc[[glue::glue("Schaefer2018_{n_parcels}Parcels_7Networks_order")]]
  
  # Format the parcellation as a \code{"xifti"}. Return it
  nv <- nrow(p$map) # 32492*2=64984
  np <- nrow(p$col) - 1
  z <- ciftiTools:::template_xifti()
  z$data$cortex_left <- p$map[seq(nv/2),,drop=FALSE]
  z$data$cortex_right <- p$map[seq(nv/2+1, nv),,drop=FALSE]
  z$meta$cifti <- list(
    intent=3007,
    brainstructures=c("left", "right"),
    names="parcels",
    labels=list()
    #misc = ...
  )
  z$meta$cifti$labels <- list(
    parcels = data.frame(
      Key = seq(0, np),
      Red = p$col$Red,
      Green = p$col$Green,
      Blue = p$col$Blue,
      Alpha = c(0, rep(1, np))
    )
  )
  rownames(z$meta$cifti$labels$parcels) <- rownames(p$col)
  z
}
