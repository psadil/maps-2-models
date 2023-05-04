
do_ptfce <- function(
    cope_files, 
    n_sub, 
    iter, 
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    storage_dir = here::here("data-raw","niis"),
    flags = "",
    resample = FALSE,
    enhance = TRUE){
  
  # bootstrap samples involve replacement!
  if (resample){
    copes <- sample(cope_files, n_sub, replace = TRUE)  
  }else{
    copes <- cope_files
  }
  
  Z <- calc_z(copes)
  
  if (enhance){
    out <- pTFCE::ptfce(
      img = Z, 
      mask = oro.nifti::readNIfTI(mask), 
      verbose = FALSE)
  }else{
    out <- list(Z = Z)
  }
  outf <- fs::path(
    storage_dir,
    glue::glue("nsub-{n_sub}_iter-{iter}_flags-{unique(flags)}"), 
    ext="qs")
  qs::qsave(out, file=outf)
  
  tibble::tibble(ptfce = outf, iter = iter, n_sub = n_sub)
}

get_active_ptfce <- function(q){
  x <- qs::qread(q)
  to_tbl0(x$Z, measure = "Z") |> 
    dplyr::filter(Z > .env$x$fwer0.05.Z) |>
    mask()
}

get_ptfce_maxes <- function(
    q,
    corrp_thresh = 0.95, 
    cluster_thresh = 0.0001, 
    mask=MNITemplate::getMNIPath("Brain_Mask", "2mm"), 
    minextent = 0){
  
  # by definition, won't see negative values when working with tfce
  # so, only need to calculate the positive clusters
  
  m <- RNifti::readNifti(mask)
  volume <- sum(m)
  x <- qs::qread(q)
  
  zfile <- fs::file_temp(ext=".nii.gz")
  RNifti::writeNifti(oro.nifti::oro2nii(x$Z), file = zfile)
  
  z <- oro.nifti::oro2nii(x$Z) |>
    neurobase::mask_img(mask=m) |>
    fslr::fsl_maths(opts=glue::glue("-thr {x$fwer0.05.Z}"))
  
  cls1 <- fslr::fslcluster(
    z,
    threshold = cluster_thresh,
    opts = glue::glue("--volume={volume} --minextent={minextent} --num=1000"))
  
  readr::read_tsv(
    cls1$olmax, col_select = c(-`...6`), show_col_types = FALSE, num_threads = 1,
    col_types = "idiii") |>
    dplyr::mutate(x=x+1,y=y+1,z=z+1) |>
    dplyr::mutate(sign = "positive") |>
    dplyr::left_join(get_sizes(cls1), by = "Cluster Index") 
}



get_ptfce_maxes_pop <- function(
    q,
    cluster_thresh = 0.0001, 
    mask=MNITemplate::getMNIPath("Brain_Mask", "2mm"), 
    minextent = 0){
  
  # by definition, won't see negative values when working with tfce
  # so, only need to calculate the positive clusters
  m <- RNifti::readNifti(mask)
  volume <- sum(m)
  x <- qs::qread(q)
  zfile <- fs::file_temp(ext = ".nii.gz")
  
  RNifti::writeNifti(oro.nifti::nii2oro(x$Z), file = zfile)
  cls1 <- fslr::fslcluster(
    zfile, 
    threshold = cluster_thresh,
    opts = glue::glue("--volume={volume} --minextent={minextent} --num=1000"))
  
  readr::read_tsv(
    cls1$olmax, col_select = c(-`...6`), show_col_types = FALSE, num_threads = 1,
    col_types = "idiii") |>
    dplyr::mutate(x=x+1,y=y+1,z=z+1) |>
    dplyr::mutate(sign = "positive") |>
    dplyr::left_join(get_sizes(cls1), by = "Cluster Index")
}



do_roi <- function(
    cope_files, 
    iter,
    at,
    n_sub = NULL,
    resample = FALSE){
  
  # bootstrap samples involve replacement!
  if (resample){
    copes <- sample(cope_files, n_sub, replace = TRUE)  
  }else{
    n_sub <- length(cope_files)
    copes <- cope_files
  }
  # zs <- stringr::str_replace(copes, "cope1.nii.gz", "zstat1.nii.gz")
  names(copes) <- seq_along(copes)

  purrr::map_dfr(copes, to_tbl, .id = "sub") |>
    dplyr::left_join(at, by = dplyr::join_by(x, y, z)) |>
    dplyr::filter(!is.na(label)) |>
    dplyr::summarise(
      Z = mean(value),
      sd = sd(value),
      .by = c(label, sub)
    ) |>
    dplyr::mutate(
      n_parcels = unique(at$n_parcels),
      iter = iter
    )
}

test_roi <- function(rois, ...,  .fwer = 0.05) {
  rois |>
    dplyr::group_nest(..., label) |>
    dplyr::mutate(
      fit = purrr::map(data, ~stats::t.test(Z ~ 1, data=.x) |> broom::tidy()),
      n_sub = purrr::map_dbl(data, nrow)
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(fit) |>
    dplyr::mutate(
      p.adjusted = stats::p.adjust(p.value, "holm"), 
      r = rank(estimate),
      .by = ...) |>
    dplyr::mutate(active = p.adjusted < .fwer) 
}
