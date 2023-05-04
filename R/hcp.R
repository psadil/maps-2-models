
not_avail <- function(){
  c(110613, 113417, 113821, 120010, 121719, 130518, 139637, 143830, 146836, 
    168139, 175035, 176239, 185038, 189652, 199958, 201515, 202820, 385046, 
    401422, 415837, 433839, 462139, 465852, 469961, 644246, 656657, 688569, 
    723141, 767464, 872764, 943862, 965367, 969476, 987983, 994273, 433839)
}

get_all_gray <- function(tfce){
  tfce <- tfce |>
    dplyr::filter(stringr::str_detect(tfce_corrp_tstat, glue::glue("flags-{ContrastName}_tfce")))
  checkmate::assert_data_frame(tfce, nrows = 1)
  
  tfce |>
    dplyr::select(Task, CopeNumber, ContrastName, n_sub, tstat, iter) |>
    dplyr::mutate(
      data2 = purrr::map2(
        tstat, n_sub, 
        ~get_pairs(stringr::str_replace(.x, "/_", "_"), .y))) |>
    tidyr::unnest(data2) |>
    dplyr::select(-tstat)
}

get_active <- function(tfce, threshold = 0.999){
  tfce <- tfce |>
    dplyr::filter(stringr::str_detect(tfce_corrp_tstat, glue::glue("flags-{ContrastName}_tfce")))
  checkmate::assert_data_frame(tfce, nrows = 1)
  gray <- to_tbl(MNITemplate::getMNISegPath(res="2mm")) |>
    dplyr::filter(value==2) |>
    dplyr::select(-value)
  
  to_tbl(tfce$tfce_corrp_tstat, measure = "p") |>
    mask() |>
    dplyr::semi_join(gray, by=c("x","y","z")) |>
    dplyr::filter(p > .env$threshold) |>
    dplyr::mutate(
      n_sub = tfce$n_sub, 
      iter = tfce$iter, 
      Task = tfce$Task, 
      ContrastName = tfce$ContrastName,
      CopeNumber = tfce$CopeNumber)
}

cor_pairwise <- function(tfce, ContrastName, n_sub, method="spearman"){
  tfce <- tfce |>
    dplyr::filter(
      .data$n_sub == .env$n_sub, 
      .data$ContrastName == .env$ContrastName) |>
    dplyr::filter(stringr::str_detect(tfce_corrp_tstat, glue::glue("flags-{ContrastName}_tfce")))
  
  tmp <- tfce |>
    dplyr::select(Task, CopeNumber, ContrastName, n_sub, tstat, iter) |>
    dplyr::mutate(
      data2 = purrr::map2(
        tstat, n_sub, 
        ~get_pairs(stringr::str_replace(.x, "/_", "_"), .y))) |>
    tidyr::unnest(data2) |>
    dplyr::select(-tstat, -sigma) |>
    tidyr::pivot_wider(names_from = iter, values_from = cope)
  
  tmp |>
    dplyr::distinct(Task, CopeNumber, ContrastName, n_sub) |>
    dplyr::mutate(
      rhos = list(
        tmp |>
          dplyr::select(tidyselect::matches("[[:digit:]]+")) |>
          corrr::correlate(method = .env$method) |>
          corrr::stretch()
      ),
      method = .env$method) |>
    tidyr::unnest(rhos)
}

cor_w_pop <- function(tfce, tfce_pop, method = "spearman"){
  tfce <- tfce |>
    dplyr::filter(stringr::str_detect(tfce_corrp_tstat, glue::glue("flags-{ContrastName}_tfce")))
  tfce_pop <- tfce_pop |>
    dplyr::semi_join(tfce, by = c("ContrastName"))
  checkmate::assert_data_frame(tfce, nrows = 1)
  checkmate::assert_data_frame(tfce_pop, nrows = 1)
  
  gray <- to_tbl(MNITemplate::getMNISegPath(res="2mm")) |>
    dplyr::filter(value==2) |>
    dplyr::select(-value)
  
  study <- get_pairs(stringr::str_replace(tfce$tstat, "/_", "_"), tfce$n_sub) |>
    mask() |>
    dplyr::semi_join(gray, by=c("x","y","z")) |>
    dplyr::select(x, y, z, study = cope)
  
  test <- get_pairs(stringr::str_replace(tfce_pop$tstat, "/_", "_"), tfce_pop$n_sub) |>
    mask() |>
    dplyr::semi_join(gray, by=c("x","y","z")) |>
    dplyr::select(x, y, z, test = cope)
  
  tfce |>
    dplyr::select(Task, CopeNumber, ContrastName, iter, n_sub) |>
    dplyr::bind_cols(
      dplyr::left_join(study, test, by = c("x", "y", "z")) |>
        dplyr::summarise(rho = cor(study, test, method = .env$method))
    ) |>
    dplyr::mutate(method = .env$method)
}

cor_w_pop2 <- function(tfce, tfce_pop, method = "spearman"){
  tfce <- tfce |>
    dplyr::filter(stringr::str_detect(tfce_corrp_tstat, glue::glue("flags-{ContrastName}_tfce")))
  tfce_pop <- tfce_pop |>
    dplyr::semi_join(tfce, by = c("ContrastName"))
  checkmate::assert_data_frame(tfce, nrows = 1)
  checkmate::assert_data_frame(tfce_pop, nrows = 1)
  
  gray <- to_tbl(MNITemplate::getMNISegPath(res="2mm")) |>
    dplyr::filter(value==2) |>
    dplyr::select(-value)
  
  study <- get_pairs(stringr::str_replace(tfce$tstat, "/_", "_"), tfce$n_sub) |>
    mask() |>
    dplyr::semi_join(gray, by=c("x","y","z")) |>
    dplyr::mutate(study = cope / sigma * correct_d(tfce$n_sub)) |>
    dplyr::select(x, y, z, study)
  
  test <- get_pairs(stringr::str_replace(tfce_pop$tstat, "/_", "_"), tfce_pop$n_sub) |>
    mask() |>
    dplyr::semi_join(gray, by=c("x","y","z")) |>
    dplyr::mutate(test = cope / sigma * correct_d(tfce_pop$n_sub)) |>
    dplyr::select(x, y, z, test)
  
  tfce |>
    dplyr::select(Task, CopeNumber, ContrastName, iter, n_sub) |>
    dplyr::bind_cols(
      dplyr::left_join(study, test, by = c("x", "y", "z")) |>
        dplyr::summarise(rho = cor(study, test, method = .env$method, use = "complete.obs"))
    ) |>
    dplyr::mutate(method = .env$method)
}

cor_w_pop_by_region <- function(tfce, tfce_pop, at, method = "spearman"){
  tfce <- tfce |>
    dplyr::filter(stringr::str_detect(tfce_corrp_tstat, glue::glue("flags-{ContrastName}_tfce")))
  tfce_pop <- tfce_pop |>
    dplyr::semi_join(tfce, by = c("ContrastName"))
  checkmate::assert_data_frame(tfce, nrows = 1)
  checkmate::assert_data_frame(tfce_pop, nrows = 1)
  
  gray <- to_tbl(MNITemplate::getMNISegPath(res="2mm")) |>
    dplyr::filter(value==2) |>
    dplyr::select(-value)
  
  study <- get_pairs(stringr::str_replace(tfce$tstat, "/_", "_"), tfce$n_sub) |>
    mask() |>
    dplyr::semi_join(gray, by=c("x","y","z")) |>
    dplyr::select(x, y, z, study = cope)
  
  test <- get_pairs(stringr::str_replace(tfce_pop$tstat, "/_", "_"), tfce_pop$n_sub) |>
    mask() |>
    dplyr::semi_join(gray, by=c("x","y","z")) |>
    dplyr::select(x, y, z, test = cope)
  
  dplyr::left_join(study, test, by = c("x", "y", "z")) |>
    dplyr::left_join(at, by = c("x", "y", "z")) |>
    dplyr::filter(!is.na(label)) |>
    dplyr::group_by(label, `Label Name`, `Network Name`, `Full component name`, hemi) |>
    dplyr::summarise(rho = cor(study, test, method = .env$method), .groups = "drop") |>
    dplyr::bind_cols(
      tfce |>
        dplyr::select(Task, CopeNumber, ContrastName, iter, n_sub)) |>
    dplyr::mutate(method = .env$method)
}


get_active0 <- function(tfce, threshold = 0.999){
  tfce <- tfce |>
    dplyr::filter(stringr::str_detect(tfce_corrp_tstat, glue::glue("flags-{ContrastName}_tfce")))
  checkmate::assert_data_frame(tfce, nrows = 1)
  gray <- to_tbl(MNITemplate::getMNISegPath(res="2mm")) |>
    dplyr::filter(value==2) |>
    dplyr::select(-value)
  
  to_tbl(stringr::str_replace(tfce$tstat, "/_", "_"), measure = "tstat") |>
    mask() |>
    dplyr::semi_join(gray, by=c("x","y","z")) |>
    dplyr::mutate(p = 1 - pt(abs(tstat), df = tfce$n_sub - 1, lower.tail = FALSE)*2) |>
    dplyr::filter(p > .env$threshold) |>
    dplyr::mutate(
      n_sub = tfce$n_sub, 
      iter = tfce$iter, 
      Task = tfce$Task, 
      ContrastName = tfce$ContrastName,
      CopeNumber = tfce$CopeNumber)
}



do_tfce2 <- function(
    cope_files, 
    n_sub, 
    iter, 
    output, 
    n = 1000, 
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    storage_dir = here::here("data-raw","niis"),
    flags = ""){
  
  # bootstrap samples involve replacement!
  copes <- sample(cope_files, n_sub, replace = TRUE)
  merged <- fsl_merge_long(copes)
  
  filestem <- fs::path(storage_dir, glue::glue("nsub-{n_sub}_iter-{iter}_flags-{flags}"))
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


do_tfce_pop2 <- function(
    cope_files, 
    n_sub, 
    iter, 
    output, 
    mask = MNITemplate::getMNIPath("Brain_Mask", "2mm"),
    storage_dir = here::here("data-raw","niis"),
    flags = ""){
  
  merged <- fsl_merge_long(cope_files)
  
  filestem <- fs::path(storage_dir, glue::glue("nsub-{n_sub}_iter-{iter}_flags-{flags}"))
  out <- fsl_randomise(merged, filestem, "-1", "-m", mask, "-R", "--uncorrp", "--glm_output")
  tibble::tibble(
    tstat = paste0(filestem, "_tstat1.nii.gz"),
    iter = iter,
    n_sub = n_sub
  )
}

get_pairs <- function(tstat, n_sub){
  cope <- to_tbl(stringr::str_replace(tstat, "_tstat1", "_glm_cope"), measure = "cope")
  varcope <- to_tbl(stringr::str_replace(tstat, "_tstat1", "_glm_varcope"), measure = "varcope") |>
    dplyr::mutate(sigma = sqrt(varcope) * sqrt(n_sub)) |>
    dplyr::select(-varcope)
  
  gray <- to_tbl(MNITemplate::getMNISegPath(res = "2mm")) |>
    dplyr::filter(value==2) |>
    dplyr::select(-value)
  
  dplyr::left_join(cope, varcope, by=c("x","y","z")) |>
    mask() |>
    dplyr::semi_join(gray, by=c("x","y","z"))
}


get_pop_d <- function(tfce_pop){
  tfce_pop |>
    dplyr::select(tstat, Task, CopeNumber, ContrastName, n_sub) |>
    dplyr::mutate(data = purrr::map2(tstat, n_sub, get_pairs)) |>
    dplyr::select(-tstat) |>
    tidyr::unnest(data) |>
    dplyr::mutate(
      d = cut(
        abs(cope / sigma), 
        right = FALSE,
        ordered = TRUE,
        breaks = c(0, .2, .5, .8, Inf),
        labels = c("0", "small", "medium", "large")))
}

format_arrow_table <- function() {
  targets::tar_format(
    read = function(path) {arrow::read_parquet(path, as_data_frame = FALSE)},
    write = function(object, path) {arrow::write_parquet(object, path, version = "2.6")},
    marshal = function(object) as.data.frame(object),
    unmarshal = function(object) arrow::Table$create(object)
  )}
