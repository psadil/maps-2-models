
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
  out <- fsl_randomise(merged, filestem, "-1", "-T", "-n", n, "-m", mask, "-R", "--uncorrp", "--glm_output")
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
    output,
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    storage_dir = here::here("data-raw", "niis"),
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


