make_atlas <- function(labelxml, nii){
  labels <- xml2::read_xml(labelxml) |> 
    xml2::as_list() |> 
    magrittr::extract2("atlas") |> 
    magrittr::extract2("data") |> 
    tibble::as_tibble(.name_repair = "unique") |> 
    tidyr::pivot_longer(cols = tidyselect::everything()) |> 
    dplyr::mutate(
      index = stringr::str_extract(name, "[[:digit:]]+") |> as.integer(), 
      label = purrr::map_chr(value, ~.x)) |>
    dplyr::select(index, label)
  
  neurobase::readnii(nii) |>
    neurobase::img_indices(add_values = TRUE) |>
    tibble::as_tibble() |>
    dplyr::rename(index = value) |>
    dplyr::filter(index > 0) |>
    dplyr::left_join(labels, by = "index") 
}


make_atlas_full <- function(){
  cortical <- make_atlas(
    fs::path(fslr::fsldir(), "data/atlases/HarvardOxford-Cortical.xml"),
    fs::path(fslr::fsldir(), "data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-2mm.nii.gz"))
  
  subcortical <- make_atlas(
    fs::path(fslr::fsldir(), "data/atlases/HarvardOxford-Subcortical.xml"),
    fs::path(fslr::fsldir(), "data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr0-2mm.nii.gz"))

  dplyr::full_join(cortical, subcortical, by = c("x", "y", "z"), suffix = c(".cortex",".subcortex")) |>
    dplyr::mutate(
      label = dplyr::case_when(
        is.na(label.cortex) & !is.na(label.subcortex) ~ label.subcortex,
        is.na(label.subcortex) & !is.na(label.cortex) ~ label.cortex,
        stringr::str_detect(label.subcortex, "Cerebral White") & !is.na(label.cortex)  ~ label.cortex,
        stringr::str_detect(label.subcortex, "Cerebral Cortex") & !is.na(label.cortex) ~ label.cortex,
        stringr::str_detect(label.subcortex, "Brain-Stem") & !is.na(label.cortex) ~ label.cortex,
        !is.na(label.subcortex) & !is.na(label.cortex) ~ label.subcortex,
        TRUE ~ NA_character_),
      hemi = dplyr::if_else(x > 45, "Left", "Right")) |>
    # if_else(is.na(label.subcortex), label.cortex, label.subcortex)) |>
    dplyr::select(-label.cortex, -label.subcortex, -starts_with("index")) |>
    dplyr::group_by(label, hemi) |>
    dplyr::mutate(
      n_voxels = dplyr::n(),
      volume = n_voxels*2^3) |>
    dplyr::ungroup() 
  
}


get_sizes <- function(cls){
  sizes <- neurobase::img_indices(cls$osize, add_values = TRUE) |>
    tibble::as_tibble() |>
    dplyr::rename(cluster_size = value)
  neurobase::img_indices(cls$oindex, add_values = TRUE) |>
    tibble::as_tibble() |>
    dplyr::filter(value > 0) |>
    dplyr::rename(`Cluster Index` = value) |>
    dplyr::left_join(sizes, by = c("x","y","z")) |>
    dplyr::distinct(`Cluster Index`, cluster_size)
  
}

get_maxes <- function(
    niifile, 
    threshold=2, 
    mask=fslr::mni_fname(mm="2", brain = TRUE, mask = TRUE), 
    at=make_atlas_full(),
    minextent = 0){
  m <- neurobase::fast_readnii(mask)
  volume <- sum(m)
  z <- neurobase::fast_readnii(niifile) |>
    neurobase::mask_img(mask=m) 
  
  cls1 <- fslr::fslcluster(
    z, 
    threshold=threshold,
    opts = glue::glue("--volume={volume} --minextent={minextent} --num=1000"))
  
  cls2 <- fslr::fslcluster(
    z * -1, 
    threshold=threshold,
    opts = glue::glue("--volume={volume} --minextent={minextent} --num=1000"))
  
  pos <- readr::read_tsv(
    cls1$olmax, col_select = c(-`...6`), show_col_types = FALSE, num_threads = 1) |>
    dplyr::mutate(x=x+1, y=y+1, z=z+1) |>
    dplyr::mutate(sign = "positive") |>
    dplyr::left_join(get_sizes(cls1), by = "Cluster Index")
  neg <- readr::read_tsv(
    cls2$olmax, col_select = c(-`...6`), show_col_types = FALSE, num_threads = 1) |>
    dplyr::mutate(x=x+1, y=y+1, z=z+1) |>
    dplyr::mutate(sign = "negative") |>
    dplyr::left_join(get_sizes(cls2), by = "Cluster Index")
  
  dplyr::bind_rows(pos, neg) |>
    dplyr::left_join(at, by = c("x","y","z"))
}

augment_distance <- function(study, reference){
  # happens when no peaks were found in study
  if(nrow(study)==0) return(reference)
  reference |>
    dplyr::filter(!is.na(.data$label)) |>
    dplyr::rowwise() |>
    dplyr::mutate(study_ind = which.min(sqrt((x - study$x)^2 + (y - study$y)^2 +(z - study$z)^2 ))) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      label.study = study$label[study_ind],
      volume.study = study$volume[study_ind],
      x.study = study$x[study_ind],
      y.study = study$y[study_ind],
      z.study = study$z[study_ind],
      d = 2*sqrt((x - x.study)^2 + (y - y.study)^2 +(z - z.study)^2 )) 
}

