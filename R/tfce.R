
correct_d <- function(N){
  # https://doi.org/10.1101/865881
  # Han Bossier1∗, Thomas E. Nichols2 & Beatrijs Moerkerke1
  exp((lgamma((N-1)/2)) - log(sqrt((N-1)/2)) - lgamma((N-2)/2))
}

d_var <- function(N, d){
  h <- correct_d(N)
  ((N - 1) * (1 + N * d^2) / (N * (N - 3)) - d^2 / h^2) * h^2
}

fsl_randomise <- function(input, output, ...){
  sys::exec_wait(
    "randomise",
    args = c("-i", input,"-o", output, ...)
  )
  fs::dir_ls(fs::path_dir(output), regexp = paste0(fs::path_file(output), "_.*(tfce|tstat1).*"))
}

fsl_merge_long <- function(infiles, output=tempfile(fileext = ".nii.gz")){
  sys::exec_wait(
    "fslmerge",
    args = c("-t", output, infiles)
  )
  output
}

do_tfcebias <- function(
    cope_files, 
    n_sub, 
    iter, 
    output, 
    n = 1000, 
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    storage_dir = here::here("data-raw","niis")){
  
  # bootstrap samples involve replacement!
  copes <- cope_files
  merged <- fsl_merge_long(copes)
  
  filestem <- fs::path(storage_dir, glue::glue("nsub-{n_sub}_iter-{iter}"))
  out <- fsl_randomise(merged, filestem, "-1", "-T", "-n", n, "-m", mask, "-R", "--uncorrp", "--glm_output")
  tibble::tibble(
    tfce_corrp_tstat = out[stringr::str_detect(out, "tfce_corrp_tstat")],
    tfce_p_tstat = out[stringr::str_detect(out, "tfce_p_tstat")],
    tfce_tstat = out[stringr::str_detect(out, "tfce_tstat")],
    tstat = fs::path(filestem, "_tstat1.nii.gz"),
    iter = iter,
    n_sub = n_sub
  )
}


# 
do_tfce <- function(
    cope_files, 
    n_sub, 
    iter, 
    output, 
    n = 1000, 
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    storage_dir = here::here("data-raw","niis")){
  
  # bootstrap samples involve replacement!
  copes <- sample(cope_files, n_sub, replace = TRUE)
  merged <- fsl_merge_long(copes)
  
  filestem <- fs::path(storage_dir, glue::glue("nsub-{n_sub}_iter-{iter}"))
  out <- fsl_randomise(merged, filestem, "-1", "-T", "-n", n, "-m", mask, "-R", "--uncorrp", "--glm_output")
  tibble::tibble(
    tfce_corrp_tstat = out[stringr::str_detect(out, "tfce_corrp_tstat")],
    tfce_p_tstat = out[stringr::str_detect(out, "tfce_p_tstat")],
    tfce_tstat = out[stringr::str_detect(out, "tfce_tstat")],
    tstat = fs::path(filestem, "_tstat1.nii.gz"),
    iter = iter,
    n_sub = n_sub
  )
}


do_tfce_pop <- function(
    cope_files, 
    n_sub, 
    iter, 
    output, 
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    storage_dir = here::here("data-raw","niis")){
  
  merged <- fsl_merge_long(cope_files)
  
  filestem <- fs::path(storage_dir, glue::glue("nsub-{n_sub}_iter-{iter}"))
  out <- fsl_randomise(merged, filestem, "-1", "-m", mask, "-R", "--uncorrp", "--glm_output")
  tibble::tibble(
    tstat = paste0(filestem, "_tstat1.nii.gz"),
    iter = iter,
    n_sub = n_sub
  )
}


get_tfce_maxes <- function(
    corrp,
    tstat,
    corrp_thresh = 0.95, 
    cluster_thresh = 0.0001, 
    mask=MNITemplate::getMNIPath("Brain_Mask", "2mm"), 
    at=make_atlas_full(),
    minextent = 0){
  
  # by definition, won't see negative values when working with tfce
  # so, only need to calculate the positive clusters
  
  m <- neurobase::fast_readnii(mask)
  volume <- sum(m)
  z <- neurobase::fast_readnii(corrp) |>
    neurobase::mask_img(mask=m) |>
    fslr::fsl_maths(opts=glue::glue("-thr {corrp_thresh} -bin -mul {tstat}"))
  
  cls1 <- fslr::fslcluster(
    z, 
    threshold = cluster_thresh,
    opts = glue::glue("--volume={volume} --minextent={minextent} --num=1000"))
  
  readr::read_tsv(
    cls1$olmax, col_select = c(-`...6`), show_col_types = FALSE, num_threads = 1,
    col_types = "idiii") |>
    dplyr::mutate(x=x+1,y=y+1,z=z+1) |>
    dplyr::mutate(sign = "positive") |>
    dplyr::left_join(get_sizes(cls1), by = "Cluster Index") |>
    dplyr::left_join(at, by = c("x","y","z"))
}


get_tfce_maxes_pop <- function(
    tstat,
    cluster_thresh = 0.0001, 
    mask=MNITemplate::getMNIPath("Brain_Mask", "2mm"), 
    at=make_atlas_full(),
    minextent = 0){
  
  # by definition, won't see negative values when working with tfce
  # so, only need to calculate the positive clusters
  m <- neurobase::fast_readnii(mask)
  volume <- sum(m)
  cls1 <- fslr::fslcluster(
    tstat, 
    threshold = cluster_thresh,
    opts = glue::glue("--volume={volume} --minextent={minextent} --num=1000"))
  
  readr::read_tsv(
    cls1$olmax, col_select = c(-`...6`), show_col_types = FALSE, num_threads = 1,
    col_types = "idiii") |>
    dplyr::mutate(x=x+1,y=y+1,z=z+1) |>
    dplyr::mutate(sign = "positive") |>
    dplyr::left_join(get_sizes(cls1), by = "Cluster Index") |>
    dplyr::left_join(at, by = c("x","y","z"))
}


fix_names <- function(d){
  d |>
    dplyr::mutate(
      tstat = stringr::str_replace(.data$tstat, "/_tstat1.nii.gz", "_tstat1.nii.gz")
    )
}

do_cor <- function(studies, ref){
  checkmate::assert_data_frame(ref, nrows = 1)
  ref_nii <- to_tbl(ref$tstat[[1]], measure = "gold") |>
    mask() |>
    dplyr::mutate(gold = gold / sqrt(ref$n_sub[[1]]) * correct_d(ref$n_sub[[1]] ) ) 
  studies |>
    fix_names() |>
    dplyr::mutate(
      d = purrr::map2(
        .data$tstat, .data$n_sub,
        ~to_tbl0(RNifti::readNifti(.x) / sqrt(.y) * correct_d(.y)) |>
          mask()
      ),
      rho = purrr::map_dbl(
        .data$d,
        ~dplyr::left_join(.x, .env$ref_nii, by = c("x", "y", "z")) |> 
          dplyr::summarise(rho=cor(value, gold)) |> 
          magrittr::use_series(rho)
      )
    )
}


to_tbl0 <- function(value, measure="value"){
  dimnames(value) <- list(
    "x" = seq_len(dim(value)[1]),
    "y" = seq_len(dim(value)[2]),
    "z" = seq_len(dim(value)[3]))
  cubelyr::as.tbl_cube(value, met_name=measure) |> 
    tibble::as_tibble()
}

mask <- function(d, mask=MNITemplate::getMNIPath("Brain_Mask", "2mm")){
  m <- to_tbl0(RNifti::readNifti(mask)) |>
    dplyr::filter(value > 0)
  dplyr::semi_join(d, m, by = c("x","y","z"))
}


get_cor2 <- function(study, ref){
  
  gold <- to_tbl(fs::path_rel(ref$tstat[[1]], "/home/ubuntu/mnt/meta/meta"), measure="gold") |>
    mask() |>
    dplyr::mutate(gold = gold / sqrt(ref$n_sub[[1]]))
  
  study |>
    dplyr::select(iter, n_sub,d) |>
    tidyr::unnest(.data$d) |>
    mask() |>
    dplyr::group_by(n_sub, x, y, z) |>
    dplyr::summarise(m = mean(value), .groups="drop") |>
    dplyr::left_join(gold)
}

get_center <- function(study, ref){
  d <- study |>
    dplyr::select(iter, n_sub, rho, d) |>
    dplyr::mutate(nsub=factor(n_sub)) 
  
  b <- d |>
    dplyr::select(-rho) |>
    tidyr::unnest(d) |>
    mask() |>
    dplyr::filter(dplyr::between(z, 20, 70), z%%5==0) |>
    dplyr::group_by(n_sub, x, y, z) |>
    dplyr::summarise(
      m = mean(value),
      s = sd(value),
      .groups="drop") |>
    dplyr::left_join(
      to_tbl(fs::path_rel(ref$tstat[[1]], "/home/ubuntu/mnt/meta/meta"), measure="gold") |>
        dplyr::mutate(gold = gold / sqrt(ref$n_sub[[1]]))) |>
    dplyr::mutate(value = m - gold, a = abs(m))
}

to_tbl <- function(file, measure="value", volumes = NULL){
  value <- RNifti::readNifti(file, volumes = volumes)
  to_tbl0(value, measure=measure)
}

make_sigmacope_bias_tbl <- function(tfce_cor, tfce_pop){
  
  pop_var <- to_tbl(
    stringr::str_replace(tfce_pop$tstat[[1]], "_tstat1.nii.gz", "_glm_sigmasqr.nii.gz"), 
    measure = "sigma_gold") |>
    mask() |>
    dplyr::mutate(sigma_gold = sqrt(sigma_gold))
  pop_cope <- to_tbl(
    stringr::str_replace(tfce_pop$tstat[[1]], "_tstat1.nii.gz", "_glm_cope.nii.gz"), 
    measure = "cope_gold") |>
    mask() 
  
  pop <- dplyr::left_join(pop_var, pop_cope, by=c("x","y","z"))
  
  tfce_cor |>
    dplyr::select(tstat, n_sub, iter) |>
    dplyr::mutate(
      sigma2file = stringr::str_replace(.data$tstat, "_tstat1.nii.gz", "_glm_sigmasqr.nii.gz"),
      sigma2 = purrr::map(sigma2file, ~to_tbl(.x, measure = "sigma") |> mask()),
      copefile = stringr::str_replace(.data$tstat, "_tstat1.nii.gz", "_glm_cope.nii.gz"),
      cope = purrr::map(copefile, ~to_tbl(.x, measure = "cope") |> mask()),
      data = purrr::map2(sigma2, cope, dplyr::left_join, by = c("x","y","z"))) |>
    dplyr::select(data, n_sub, iter) |>
    tidyr::unnest(data) |>
    dplyr::group_by(x, y, z, n_sub) |>
    dplyr::summarise(
      sigma = mean(sigma),
      cope = mean(cope),
      .groups = "drop") |>
    dplyr::mutate(
      sigma = sqrt(sigma),
      src = "study") |>
    dplyr::left_join(pop, by = c("x","y","z"))
}


make_pop_d <- function(tfce_pop){
  to_tbl(tfce_pop$tstat, measure = "gold") |> 
    mask() |> 
    dplyr::mutate(
      gold = gold / sqrt(tfce_pop$n_sub[[1]]) * correct_d(tfce_pop$n_sub[[1]]),
      gold_cut = cut(
        gold, 
        include.lowest = TRUE,
        breaks=quantile(gold, seq(0, 1, by=.1)),
        glue::glue(
          "({round(quantile(gold, seq(0, 1, by=.1))[1:10],2)}, {round(quantile(gold, seq(0, 1, by=.1))[2:11], 2)}]")))  
}


make_iters <- function(tfce_cor, pop_d){
  
  tfce_cor |>
    dplyr::select(iter, n_sub, d) |>
    tidyr::unnest(d) |>
    dplyr::right_join(pop_d) |>
    dplyr::group_by(gold_cut, n_sub, iter) |>
    dplyr::summarise(
      s = sd(value),
      N = dplyr::n(),
      q.05 = quantile(value, 0.05),
      q.50 = median(value),
      average = mean(value),
      q.95 = quantile(value, 0.95),
      rho = cor(gold, value),
      rhos = cor(gold, value, method = "spearman"),
      gold = mean(gold),
      .groups = "drop") 
}

