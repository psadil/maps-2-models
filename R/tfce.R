
fsl_randomise <- function(input, output, ...) {
  sys::exec_wait(
    "randomise",
    args = c("-i", input, "-o", output, ...)
  )
  fs::dir_ls(fs::path_dir(output), regexp = paste0(fs::path_file(output), "_.*(tfce|tstat1).*"))
}

fsl_merge_long <- function(infiles, output = tempfile(fileext = ".nii.gz")) {
  sys::exec_wait(
    "fslmerge",
    args = c("-t", output, infiles)
  )
  output
}

do_tfce2 <- function(
    copes,
    n_sub,
    iter,
    n = 1000,
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    storage_dir = here::here("data-raw", "niis"),
    flags = "") {
  
  merged <- fsl_merge_long(copes)
  
  filestem <- fs::path(storage_dir, glue::glue("nsub-{n_sub}_iter-{iter}_flags-{flags}"))
  out <- fsl_randomise(
    merged, 
    filestem, 
    "-1", 
    "-T", 
    "-n", n, 
    "-m", mask, 
    "-R", "--uncorrp", "--glm_output")
  tibble::tibble(
    tfce_corrp_tstat = out[stringr::str_detect(out, "tfce_corrp_tstat")],
    tfce_p_tstat = out[stringr::str_detect(out, "tfce_p_tstat")],
    tfce_tstat = out[stringr::str_detect(out, "tfce_tstat")],
    tstat = fs::path(stringr::str_c(filestem, "_tstat1.nii.gz"))
  )
}

do_tfce_pop2 <- function(
    cope_files,
    n_sub,
    iter,
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    storage_dir = Sys.getenv("NIIDIR"),
    flags = "") {
  merged <- fsl_merge_long(cope_files)
  
  filestem <- fs::path(storage_dir, glue::glue("nsub-{n_sub}_iter-{iter}_flags-{flags}"))
  out <- fsl_randomise(merged, filestem, "-1", "-m", mask, "-R", "--uncorrp", "--glm_output")
  tibble::tibble(
    tstat = paste0(filestem, "_tstat1.nii.gz"),
    iter = iter,
    n_sub = n_sub
  )
}

.do_palm_volume <- function(hcp_samples, filestem, n=10000){
  copes <- hcp_samples |> 
    purrr::pluck("value") |> 
    stringr::str_c(collapse = " ")
  
  glue::glue("palm2 -f {filestem} {copes}")
}

.do_palm_surface <- function(hcp_samples, filestem, n=10000){
  
  copes <- hcp_samples |> 
    purrr::pluck("value") |> 
    stringr::str_c(collapse = " ")
  
  glue::glue("./tools/dopalm -f {filestem} {copes}")
}


get_tfce_cmd <- function(hcp_samples, test, n = 10000) {
  
  hcp_samples |>
    dplyr::select(-tar_group) |>
    dplyr::left_join(test, by = dplyr::join_by(Task, CopeNumber, sub)) |>
    dplyr::group_nest(type, n_sub, iter, Task, CopeNumber) |>
    dplyr::mutate(
      fit = purrr::pmap(
        list(
          .data=data, 
          type=type, 
          n_sub=n_sub, 
          iter=iter, 
          Task=Task, 
          CopeNumber=CopeNumber),
        .get_tfce_cmd
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(fit)
}



.get_tfce_cmd <- function(
    .data, 
    type,
    n_sub,
    iter,
    Task,
    CopeNumber,
    n = 10000) {
  
  stem <- glue::glue("nsub-{n_sub}_iter-{iter}_task-{Task}_cope-{CopeNumber}_type-{type}")
  prefix <- fs::path(Sys.getenv("NIIDIR"), stem)
  
  if (type%in%c("VOL","UKB")){
    cmd <- .data |>
      dplyr::rename(value=VOL) |>
      .do_palm_volume(n=n, filestem=prefix)
  }else{
    if (type=="MSMALL"){
      cmd <- .data |>
        dplyr::rename(value=MSMALL) |>
        .do_palm_surface(n=n, filestem=prefix)
    }else{
      cmd <- .data |>
        dplyr::rename(value=SURFACE) |>
        .do_palm_surface(n=n, filestem=prefix)
    }
  }
  
  dplyr::tibble(filestem=stem, cmd=cmd)
}

get_tfce_pop <- function(test, n = 1) {
  
  test |>
    tidyr::pivot_longer(c(MSMALL, SURFACE, VOL), names_to = "type") |>
    dplyr::distinct(Task, CopeNumber, type, sub) |>
    dplyr::mutate(n_sub = dplyr::n_distinct(sub), .by = c(Task, CopeNumber, type)) |>
    dplyr::mutate(
      group = interaction(Task, CopeNumber, type),
      tar_group=1,
      iter=0) |>
    dplyr::group_nest(group) |>
    dplyr::mutate(
      cmd = purrr::map(data, get_tfce_cmd, test=test)
    ) |>
    dplyr::select(cmd) |>
    tidyr::unnest(cmd)
}

.get_palm_maxes_vol <- function(samples, cluster_thresh = 0.0001){
  filestem <- samples$filestem
  tstat <- fs::path(Sys.getenv("NIIDIR"), glue::glue("{filestem}_tfce_tstat_c1.nii"))
  fwep <- fs::path(Sys.getenv("NIIDIR"), glue::glue("{filestem}_tfce_tstat_fwep_c1.nii"))
  
  zfile <- fs::file_temp(ext = ".nii.gz")
  mask <- fslr::fslthresh(fwep, uthresh = 0.05) |>
    fslr::fslbin()
  
  cls1 <- fslr::fslmask(tstat, mask) |>
    fslr::fslmask(MNITemplate::getMNIPath("Brain_Mask", "2mm")) |>
    fslr::fslcluster(
      threshold = cluster_thresh,
      opts = glue::glue("--minextent=0 --num=1000")
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

get_palm_maxes <- function(
    samples,
    corrp_thresh = 0.95,
    cluster_thresh = 0.0001) {
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


