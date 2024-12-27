calc_z <- function(cope_files) {
  copes <- RNifti::readNifti(cope_files) |> simplify2array()
  stopifnot(length(dim(copes)) == 4)
  n <- dim(copes)[4]
  means <- apply(copes, 1:3, mean)
  sds <- apply(copes, 1:3, sd)
  t_stat <- sqrt(n) * means / sds
  z_stat <- apply(
    t_stat, 1:3,
    FUN = function(x) {
      pt(x, n - 1, log.p = TRUE, lower.tail = FALSE) |>
        qnorm(log.p = TRUE, lower.tail = FALSE)
    }
  )
  z_stat[is.na(z_stat)] <- 0
  
  neurobase::niftiarr(neurobase::readnii(cope_files[[1]]), z_stat)
}

do_ptfce <- function(
    cope_files,
    n_sub,
    iter,
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    storage_dir = here::here("data-raw", "niis"),
    flags = "",
    resample = FALSE,
    enhance = TRUE) {
  # bootstrap samples involve replacement!
  if (resample) {
    copes <- sample(cope_files, n_sub, replace = TRUE)
  } else {
    copes <- cope_files
  }
  
  Z <- calc_z(copes)
  
  if (enhance) {
    out <- pTFCE::ptfce(
      img = Z,
      mask = oro.nifti::readNIfTI(mask),
      verbose = FALSE
    )
  } else {
    out <- list()
  }
  out$Z_raw <- Z
  outf <- fs::path(
    storage_dir,
    glue::glue("nsub-{n_sub}_iter-{iter}_flags-{unique(flags)}"),
    ext = "qs"
  )
  qs::qsave(out, file = outf)
  
  tibble::tibble(ptfce = outf, iter = iter, n_sub = n_sub, copes = list(copes))
}


get_active_ptfce <- function(q) {
  x <- qs::qread(q)
  to_tbl0(x$Z, measure = "Z") |>
    dplyr::filter(Z > .env$x$fwer0.05.Z) |>
    mask()
}

get_ptfce_maxes <- function(
    row,
    do_fwe_correction = TRUE,
    cluster_thresh = 0.0001,
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    minextent = 0) {
  # by definition, won't see negative values when working with tfce
  # so, only need to calculate the positive clusters
  
  if (is.na(row$ptfce)){
    return(tibble::tibble())
  }
  
  m <- RNifti::readNifti(mask)
  volume <- sum(m)
  
  q <- fs::path(fs::path(Sys.getenv("NIIDIR"), row$ptfce[[1]], ext = "qs"))
  x <- qs::qread(q)
  
  z <- oro.nifti::oro2nii(x$Z) |>
    neurobase::mask_img(mask = m) 

  if (do_fwe_correction){
    z <- z |>
      fslr::fsl_thresh(thresh= x$fwer0.05.Z)
  }
  
  cls1 <- fslr::fslcluster(
    z,
    threshold = cluster_thresh,
    opts = glue::glue("--volume={volume} --minextent={minextent} --num=1000")
  )
  
  readr::read_tsv(
    cls1$olmax,
    col_select = c(-`...6`), show_col_types = FALSE, num_threads = 1,
    col_types = "idiii"
  ) |>
    dplyr::mutate(x = x + 1, y = y + 1, z = z + 1) |>
    dplyr::mutate(sign = "positive") |>
    dplyr::left_join(get_sizes(cls1), by = "Cluster Index") |>
    dplyr::mutate(
      Task = unique(row$Task),
      iter = unique(row$iter),
      n_sub = unique(row$n_sub),
      CopeNumber = unique(row$CopeNumber),
      fwe_correction=do_fwe_correction
    )
}



get_ptfce_maxes_pop <- function(
    row,
    cluster_thresh = 0.0001,
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    minextent = 0) {
  
  # by definition, won't see negative values when working with tfce
  # so, only need to calculate the positive clusters
  m <- RNifti::readNifti(mask)
  volume <- sum(m)
  
  q <- fs::path(fs::path(Sys.getenv("NIIDIR"), row$ptfce[[1]], ext = "qs"))
  x <- qs::qread(q)
  zfile <- fs::file_temp(ext = ".nii.gz")
  
  RNifti::writeNifti(oro.nifti::nii2oro(x$Z), file = zfile)
  cls1 <- fslr::fslcluster(
    zfile,
    threshold = cluster_thresh,
    opts = glue::glue("--volume={volume} --minextent={minextent} --num=1000")
  )
  
  readr::read_tsv(
    cls1$olmax,
    col_select = c(-`...6`), show_col_types = FALSE, num_threads = 1,
    col_types = "idiii"
  ) |>
    dplyr::mutate(x = x + 1, y = y + 1, z = z + 1) |>
    dplyr::mutate(sign = "positive") |>
    dplyr::left_join(get_sizes(cls1), by = "Cluster Index") |>
    dplyr::mutate(
      Task = unique(row$Task),
      iter = unique(row$iter),
      n_sub = unique(row$n_sub),
      CopeNumber = unique(row$CopeNumber)
    )
}

get_ptfce_maxes_sub <- function(
    zstat,
    cluster_thresh = 0.0001,
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    minextent = 0) {
  m <- RNifti::readNifti(mask)
  volume <- sum(m)
  
  cls1 <- fslr::fslcluster(
    zstat,
    threshold = cluster_thresh,
    opts = glue::glue("--volume={volume} --minextent={minextent} --num=5000")
  )
  
  readr::read_tsv(
    cls1$olmax,
    col_select = c(-`...6`), show_col_types = FALSE, num_threads = 1,
    col_types = "idiii"
  ) |>
    dplyr::mutate(x = x + 1, y = y + 1, z = z + 1) |>
    dplyr::mutate(sign = "positive") |>
    dplyr::left_join(get_sizes(cls1), by = "Cluster Index")
}




do_ptfce2 <- function(
    hcp_samples,
    test,
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    enhance = TRUE) {
  
  type <- unique(hcp_samples$type)
  if (!type == "VOL"){
    out <- dplyr::distinct(hcp_samples, Task, CopeNumber, n_sub, iter) |>
      dplyr::mutate(ptfce = NA_character_)
    return(out)
  }
  
  copes <- test |>
    dplyr::select(-SURFACE, -MSMALL, -tar_group) |>
    dplyr::semi_join(hcp_samples, by = dplyr::join_by(Task, sub, CopeNumber))
  
  Z <- calc_z(copes$VOL)
  
  if (enhance) {
    out <- pTFCE::ptfce(
      img = Z,
      mask = oro.nifti::readNIfTI(mask),
      verbose = FALSE
    )
  } else {
    out <- list()
  }
  out$Z_raw <- Z
  
  n_sub <- unique(hcp_samples$n_sub)
  iter <- unique(hcp_samples$iter)
  Task <- unique(hcp_samples$Task)
  CopeNumber <- unique(hcp_samples$CopeNumber)
  
  stem <- glue::glue("nsub-{n_sub}_iter-{iter}_task-{Task}_cope-{CopeNumber}_type-{type}")
  
  outf <- fs::path(Sys.getenv("NIIDIR"), stem, ext = "qs")
  qs::qsave(out, file = outf)
  
  dplyr::distinct(hcp_samples, Task, CopeNumber, n_sub, iter) |>
    dplyr::mutate(ptfce = as.character(stem))
}

do_ptfce_pop <- function(test){
  dplyr::distinct(test, Task, CopeNumber, sub, tar_group) |>
    dplyr::mutate(iter = 0, type = "VOL") |>
    dplyr::mutate(n_sub = dplyr::n(), .by = tar_group) |>
    dplyr::group_nest(tar_group) |>
    dplyr::mutate(
      out = purrr::map(
        data,
        ~do_ptfce2(hcp_samples = .x, test=test, enhance = FALSE)
      )
    ) |>
    dplyr::select(-tar_group, -data) |>
    tidyr::unnest(out)
}
