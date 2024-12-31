
read_cifti_with_mwall <- function(file){
  ciftiTools::read_cifti(file) |>
    ciftiTools::move_from_mwall()
}

create_dense_from_template <- function(template, volume_all, cortex_left, cortex_right, out=NULL){
  if (is.null(out)){
    out <- tempfile(fileext = ".dscalar.nii")
  }
  cmd <- glue::glue("-cifti-create-dense-from-template {template} {out} -volume-all {volume_all} -metric CORTEX_LEFT {cortex_left} -metric CORTEX_RIGHT {cortex_right}")
  check <- ciftiTools::run_wb_cmd(cmd, intern = FALSE)
  if (!check){
    stop("Failed to create dense from template. cmd: ", cmd)
  }
  out
}

mask_cifti_by_logp <- function(img, fwe, out=NULL, threshold=-log(0.05)){
  if (is.null(out)){
    out <- tempfile(fileext = ".dscalar.nii")
  }
  fwe_cii <- ciftiTools::read_cifti(fwe)
  mask <- fwe_cii > threshold
  tstat_cii <- ciftiTools::read_cifti(img)
  masked <- tstat_cii * mask
  ciftiTools::write_cifti(masked, out)
}

surface_average <- function(surfs, out=NULL){
  if (is.null(out)){
    out <- tempfile(fileext = ".surf.gii")
  }
  to_merge <- stringr::str_c("-surf ", surfs) |> paste0(collapse = " ")
  check <- ciftiTools::run_wb_cmd(
    glue::glue("-surface-average {out} {to_merge}"),
    intern = FALSE
  )
  if (!check){
    stop("Failed to average surface")
  }
  
  out
}

metric_merge <- function(metric, out = NULL){
  if (is.null(out)){
    out <- tempfile(fileext = ".func.gii")
  }
  to_merge <- stringr::str_c("-metric ", metric) |> paste0(collapse = " ")
  check <- ciftiTools::run_wb_cmd(
    glue::glue("-metric-merge {out} {to_merge}"),
    intern = FALSE
  )
  if (!check){
    stop("Failed to merge metrics")
  }
  
  out
}

metric_reduce <- function(metric_in, metric_out=NULL, .f="MEAN"){
  if (is.null(metric_out)){
    metric_out <- tempfile(fileext = ".shape.gii")
  }
  check <- ciftiTools::run_wb_cmd(
    glue::glue("-metric-reduce {metric_in} {.f} {metric_out}"),
    intern = FALSE
  )
  if (!check){
    stop("Failed to reduce metrics")
  }
  
  metric_out
}

get_maxima <- function(
    cifti, 
    left_surface, 
    right_surface, 
    cifti_out = NULL, 
    surface_distance=1,
    volume_distance=1, 
    direction="COLUMN"
){
  if (is.null(cifti_out)){
    cifti_out <- tempfile(fileext=".dscalar.nii")
  }
  check <- ciftiTools::run_wb_cmd(
    glue::glue("-cifti-extrema {cifti} {surface_distance} {volume_distance} {direction} {cifti_out} -left-surface {left_surface} -right-surface {right_surface} -only-maxima"),
    intern = FALSE
  )
  if (!check){
    stop("Failed to find extrema")
  }
  
  cifti_out
}

get_gold_cohens <- function(nii_gold, out=NULL){
  if (is.null(out)){
    out <- tempfile(fileext = ".dscalar.nii")
  }
  # get gold standard extrema
  xii_gold <- ciftiTools::read_cifti(nii_gold)
  
  nii_gold_avg <- ciftiTools::apply_xifti(xii_gold, margin=1, mean)
  nii_gold_sd <- ciftiTools::apply_xifti(xii_gold, margin=1, sd)
  
  nii_cohens <- nii_gold_avg / nii_gold_sd
  
  ciftiTools::write_cifti(nii_cohens, out)
}


#' get distances from gold
get_1_min_distance <- function(index, extrema, surface, corrected_areas, structure){
  distances <- tempfile(fileext = ".shape.gii")
  
  # -1 because wb_command indexing starts at 0
  # and all stored values are starting at 1
  check <- ciftiTools::run_wb_cmd(
    glue::glue("-surface-geodesic-distance {surface} {index-1} {distances} -corrected-areas {corrected_areas}"),
    intern = FALSE
  )
  if (!check){
    stop("Failed to calculate geodesic distance")
  }
  
  
  maxes <- extrema |>
    dplyr::filter(stringr::str_detect(.env$structure, .data$structure)) |>
    magrittr::use_series("index")
  
  l_distances <- freesurferformats::read.fs.morph(distances)
  
  tibble::tibble(d=l_distances[maxes], index.study=maxes) |>
    dplyr::slice_min(order_by = d, n=1, with_ties = FALSE)
}

get_cifti_vox_peaks <- function(file, extra=NULL){
  xii <- read_cifti_with_mwall(file)
  
  # https://github.com/mandymejia/ciftiTools?tab=readme-ov-file#how-do-i-get-voxelindicesijk-or-the-mni-coordinates-for-the-subcortex
  # do not subtracting 1 because these indexes will be used in R
  VoxIJK <- which(xii$meta$subcort$mask, arr.ind=TRUE)
  
  peaks_xyz <- VoxIJK[xii$data$subcort==1, 1:3] |>
    matrix(ncol=3)
  colnames(peaks_xyz) <- c("x", "y", "z")
  
  if(!is.null(extra)){
    extra_xii <- read_cifti_with_mwall(extra)
    value <- extra_xii$data$subcort[xii$data$subcort==1,]
  }else{
    value <- NA
  }
  
  tibble::as_tibble(peaks_xyz) |>
    dplyr::mutate(value=value)
  
}

get_min_distance_vol <- function(src, extrema, extrema_gold, vox_mm = 2){

  peaks_gold <- get_cifti_vox_peaks(extrema_gold, extra=src) 
  augment_distance(study = extrema, reference = peaks_gold, vox_mm = vox_mm) |>
    dplyr::mutate(structure="subcort") |>
    dplyr::select(-study_ind)
}

get_min_distance_surf <- function(
    extrema, 
    extrema_gold, 
    mid_merged, 
    va_mean, 
    structure, 
    cohens_gold
){
  
  ee_gold <- read_cifti_with_mwall(extrema_gold)
  cohens_gold_xii <- read_cifti_with_mwall(cohens_gold)
  
  tibble::tibble(
    index = which(ee_gold$data[[structure]]==1),
  ) |>
    dplyr::mutate(
      distances = purrr::map(
        index,
        get_1_min_distance,
        extrema=extrema,
        surface=mid_merged,
        corrected_areas=va_mean,
        structure=structure
      ),
      value = cohens_gold_xii$data[[structure]][index],
      structure=structure
    ) |>
    tidyr::unnest(distances)
}

get_masked_tstat <- function(filestem){
  #' threshold tfce_tstat_c1 by tfce_tstat_fwep_c1 (neg)
  
  niidir <- Sys.getenv("NIIDIR")
  palmdir <- Sys.getenv("PALMDIR")
  template <- glue::glue("{palmdir}/{filestem}.dtseries.nii")

  fwe <- create_dense_from_template(
    template = template, 
    volume_all = glue::glue("{niidir}/{filestem}_tfce_tstat_fwep_c1.nii"), 
    cortex_left = glue::glue("{niidir}/{filestem}_L_tfce_tstat_fwep_c1.gii"), 
    cortex_right = glue::glue("{niidir}/{filestem}_R_tfce_tstat_fwep_c1.gii")
  )
  
  tstat <- create_dense_from_template(
    template = template,
    volume_all = glue::glue("{niidir}/{filestem}_tfce_tstat_c1.nii"), 
    cortex_left = glue::glue("{niidir}/{filestem}_L_tfce_tstat_c1.gii"),
    cortex_right = glue::glue("{niidir}/{filestem}_R_tfce_tstat_c1.gii"),
  )
  
  mask_cifti_by_logp(img = tstat, fwe = fwe)
}


get_cifti_augmented <- function(extrema, filestem, filestem_gold){
  
  niidir <- Sys.getenv("NIIDIR")
  palmdir <- Sys.getenv("PALMDIR")
  
  # always check distances with gold surfaces and areas
  l_mid_gold <- glue::glue("{palmdir}/{filestem_gold}_L_midthickness.32k_fs_LR.gii")
  r_mid_gold <- glue::glue("{palmdir}/{filestem_gold}_R_midthickness.32k_fs_LR.gii")

  nii_gold <- glue::glue("{palmdir}/{filestem_gold}.dtseries.nii")
  cohens_gold <- get_gold_cohens(nii_gold)
  extrema_gold <- get_maxima(
    cifti=cohens_gold,
    left_surface = l_mid_gold,
    right_surface = r_mid_gold
  )
  
  left_distances <- get_min_distance_surf(
    extrema=extrema, 
    extrema_gold=extrema_gold, 
    mid_merged=l_mid_gold, 
    va_mean=glue::glue("{palmdir}/{filestem_gold}_L_area.32k_fs_LR.gii"), 
    structure="cortex_left",
    cohens_gold=cohens_gold
  )
  right_distances <- get_min_distance_surf(
    extrema=extrema,
    extrema_gold=extrema_gold,
    mid_merged=r_mid_gold, 
    va_mean=glue::glue("{palmdir}/{filestem_gold}_R_area.32k_fs_LR.gii"),
    structure="cortex_right",
    cohens_gold=cohens_gold
  )
  
  vol_distances <- get_min_distance_vol(
    cohens_gold,
    extrema=dplyr::filter(extrema, structure=="subcort"),
    extrema_gold=extrema_gold
  )
  
  dplyr::bind_rows(left_distances, right_distances, vol_distances)
}

get_cifti_augmented2 <- function(extrema, tfce_pop){
  
  extrema |>
    dplyr::select(-value) |>
    dplyr::group_nest(type, Task, CopeNumber, n_sub, filestem) |>
    dplyr::left_join(
      dplyr::select(tfce_pop, -n_sub, -iter, -cmd), 
      by = dplyr::join_by(type, Task, CopeNumber)) |>
    dplyr::mutate(
      distances = purrr::pmap(
        list(extrema=data, filestem=filestem.x, filestem_gold=filestem.y),
        get_cifti_augmented
      )
    ) |>
    dplyr::select(-data, -tidyselect::starts_with("filestem")) |>
    tidyr::unnest(distances)
}


masked_cifti_to_tbl <- function(xii, extrema, structure){
  tibble::tibble(
    index = which(extrema$data[[structure]]==1),
  ) |>
    dplyr::mutate(value=xii$data[[structure]][index], structure=structure)
}


get_cifti_peaks <- function(filestem){
  
  niidir <- Sys.getenv("NIIDIR")
  palmdir <- Sys.getenv("PALMDIR")
  
  cifti <- get_masked_tstat(filestem)
  
  maxima_file <- get_maxima(
    cifti = cifti,
    left_surface = glue::glue("{palmdir}/{filestem}_L_midthickness.32k_fs_LR.gii"), 
    right_surface = glue::glue("{palmdir}/{filestem}_R_midthickness.32k_fs_LR.gii"),
  ) 
  maxima <- maxima_file |>
    read_cifti_with_mwall()
  
  xii <- read_cifti_with_mwall(cifti)
  
  dplyr::bind_rows(
    masked_cifti_to_tbl(xii, maxima, "cortex_left"),
    masked_cifti_to_tbl(xii, maxima, "cortex_right"),
    get_cifti_vox_peaks(maxima_file, extra = cifti) |>
      dplyr::mutate(structure="subcort")
  )
  
}

get_nifti_peaks <- function(
    filestem, 
    cluster_thresh = 0.0001,
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    minextent = 0
){
  
  m <- RNifti::readNifti(mask)
  volume <- sum(m)
  
  niidir <- Sys.getenv("NIIDIR")
  tstat <- glue::glue("{niidir}/{filestem}_tfce_tstat_c1.nii")
  fwe <- glue::glue("{niidir}/{filestem}_tfce_tstat_fwep_c1.nii")
  
  fwe_nii <- RNifti::readNifti(fwe) 
  fwe_mask <- RNifti::asNifti(fwe_nii > -log(0.05), reference = fwe_nii)
  
  tstat_nii <- RNifti::readNifti(tstat) |>
    neurobase::mask_img(mask = m) |>
    neurobase::mask_img(mask = fwe_mask) 
  
  cls1 <- fslr::fslcluster(
    tstat_nii,
    threshold = cluster_thresh,
    opts = glue::glue("--volume={volume} --minextent={minextent} --num=1000")
  )
  
  readr::read_tsv(
    cls1$olmax,
    col_select = c(-`...6`), 
    show_col_types = FALSE, 
    num_threads = 1,
    col_types = "idiii"
  ) |>
    dplyr::mutate(x = x + 1, y = y + 1, z = z + 1) |>
    dplyr::left_join(get_sizes(cls1), by = "Cluster Index")
  
}

get_study_peaks_cifti <- function(tfce_row){
  type <- unique(tfce_row$type)
  
  if (type %in% c("SURFACE", "MSMALL")){
    out <- tfce_row |>
      dplyr::mutate(
        extrema = purrr::map(filestem, get_cifti_peaks) 
      )
  }else{
    out <- tfce_row |>
      dplyr::mutate(
        extrema = purrr::map(filestem, get_nifti_peaks) 
      )
  }
  
  out |>
    dplyr::select(-cmd) |>
    tidyr::unnest(extrema)
}

