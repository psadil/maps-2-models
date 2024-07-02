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

do_ptfce_sub <- function(
    zstat,
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    storage_dir = here::here("data-raw", "niis"),
    flags = "",
    enhance = TRUE) {
  Z <- oro.nifti::readNIfTI(zstat)

  if (enhance) {
    out <- pTFCE::ptfce(
      img = zstat,
      mask = oro.nifti::readNIfTI(mask),
      verbose = FALSE
    )
  } else {
    out <- list()
  }
  out$Z_raw <- Z
  outf <- fs::path(
    storage_dir,
    glue::glue("sub-{sub}_flags-{unique(flags)}"),
    ext = "qs"
  )
  qs::qsave(out, file = outf)

  tibble::tibble(ptfce = outf, iter = 0, n_sub = 1, copes = list(copes))
}

get_active_ptfce <- function(q) {
  x <- qs::qread(q)
  to_tbl0(x$Z, measure = "Z") |>
    dplyr::filter(Z > .env$x$fwer0.05.Z) |>
    mask()
}

get_ptfce_maxes <- function(
    q,
    corrp_thresh = 0.95,
    cluster_thresh = 0.0001,
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    minextent = 0) {
  # by definition, won't see negative values when working with tfce
  # so, only need to calculate the positive clusters

  m <- RNifti::readNifti(mask)
  volume <- sum(m)
  x <- qs::qread(q)

  zfile <- fs::file_temp(ext = ".nii.gz")
  RNifti::writeNifti(oro.nifti::oro2nii(x$Z), file = zfile)

  z <- oro.nifti::oro2nii(x$Z) |>
    neurobase::mask_img(mask = m) |>
    fslr::fsl_maths(opts = glue::glue("-thr {x$fwer0.05.Z}"))

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
    dplyr::left_join(get_sizes(cls1), by = "Cluster Index")
}



get_ptfce_maxes_pop <- function(
    q,
    cluster_thresh = 0.0001,
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    minextent = 0) {
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
    opts = glue::glue("--volume={volume} --minextent={minextent} --num=1000")
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


do_roi <- function(
    cope_files,
    iter,
    at,
    n_sub = NULL,
    resample = FALSE) {
  # bootstrap samples involve replacement!
  if (resample) {
    copes <- sample(cope_files, n_sub, replace = TRUE)
  } else {
    n_sub <- length(cope_files)
    copes <- cope_files
  }
  # zs <- stringr::str_replace(copes, "cope1.nii.gz", "zstat1.nii.gz")
  names(copes) <- seq_along(copes)

  purrr::map_dfr(
    copes,
    ~ to_tbl(.x) |>
      dplyr::left_join(at, by = dplyr::join_by(x, y, z)) |>
      dplyr::filter(!is.na(label)),
    .id = "sub"
  ) |>
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

test_roi <- function(rois, ..., .fwer = 0.05) {
  rois |>
    dplyr::group_nest(..., label) |>
    dplyr::mutate(
      fit = purrr::map(data, ~ stats::t.test(Z ~ 1, data = .x) |> broom::tidy()),
      n_sub = purrr::map_dbl(data, nrow)
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(fit) |>
    dplyr::mutate(
      p.adjusted = stats::p.adjust(p.value, "holm"),
      r = rank(estimate),
      .by = ...
    ) |>
    dplyr::mutate(active = p.adjusted < .fwer)
}


cor_pairwise_ptfce <- function(tfce, ContrastName, n_sub, method = "spearman") {
  tfce <- tfce |>
    dplyr::filter(
      .data$n_sub == .env$n_sub,
      .data$ContrastName == .env$ContrastName
    )

  tmp <- tfce |>
    dplyr::select(Task, CopeNumber, ContrastName, n_sub, ptfce, iter) |>
    dplyr::mutate(
      data2 = purrr::map(
        ptfce,
        ~ to_tbl0(qs::qread(.x)$Z) |> mask()
      )
    ) |>
    tidyr::unnest(data2) |>
    dplyr::select(-ptfce) |>
    tidyr::pivot_wider(names_from = iter, values_from = value)

  tmp |>
    dplyr::distinct(Task, CopeNumber, ContrastName, n_sub) |>
    dplyr::mutate(
      rhos = list(
        tmp |>
          dplyr::select(tidyselect::matches("[[:digit:]]+")) |>
          corrr::correlate(method = .env$method, quiet = TRUE) |>
          corrr::stretch()
      ),
      method = .env$method
    ) |>
    tidyr::unnest(rhos)
}
