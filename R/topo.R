.get_glm <- function(
  files,
  Task,
  CopeNumber,
  type,
  iter,
  n_sub
) {
  if (type %in% c("MSMALL", "SURFACE")) {
    xiis <- purrr::map(files$src, ciftiTools::read_cifti) |>
      ciftiTools::merge_xifti(xifti_list = _)
    avgs <- ciftiTools::apply_xifti(xiis, margin = 1, FUN = mean) |>
      as.matrix() |>
      as.vector()
    sds <- ciftiTools::apply_xifti(xiis, margin = 1, FUN = sd) |>
      as.matrix() |>
      as.vector()
  } else {
    merged <- fsl_merge_long(files$src)
    avgs <- fslr::fslmaths(merged, opts = "-Tmean") |>
      to_tbl0() |>
      mask() |>
      mask_gray() |>
      mask_atlas() |>
      magrittr::use_series("value")

    sds <- fslr::fslmaths(merged, opts = "-Tstd") |>
      to_tbl0() |>
      mask() |>
      mask_gray() |>
      mask_atlas() |>
      magrittr::use_series("value")
  }
  out_file <- fs::path(
    Sys.getenv("NIIDIR"),
    glue::glue(
      "nsub-{n_sub}_iter-{iter}_task-{Task}_cope-{CopeNumber}_type-{type}_glm.parquet"
    )
  )
  tibble::tibble(pe = avgs, sigma = sds) |>
    dplyr::mutate(
      index = 1:dplyr::n(),
      z = (sqrt(n_sub) * avgs / sds) |> # t_stat -> z
        pt(n_sub - 1, log.p = TRUE, lower.tail = FALSE) |>
        qnorm(log.p = TRUE, lower.tail = FALSE)
    ) |>
    arrow::write_parquet(out_file, compression = "zstd")

  tibble::tibble(
    Task = Task,
    CopeNumber = CopeNumber,
    type = type,
    iter = iter,
    n_sub = n_sub,
    glm = out_file
  )
}


get_glm <- function(test, samples) {
  if ("MSMALL" %in% names(test)) {
    files <- test |>
      tidyr::pivot_longer(
        c(MSMALL, SURFACE, VOL),
        names_to = "type",
        values_to = "src"
      ) |>
      dplyr::select(-tar_group)
  } else {
    files <- test |>
      tidyr::pivot_longer(UKB, names_to = "type", values_to = "src") |>
      dplyr::mutate(
        src = stringr::str_replace(
          src,
          "/fastscratch/myscratch/pssadil/ukb_mni/derivatives",
          Sys.getenv("UKBMNI")
        )
      )
  }
  files <- files |>
    na.omit() |>
    dplyr::right_join(
      samples,
      by = dplyr::join_by(Task, sub, CopeNumber, type)
    )

  Task <- unique(samples$Task)
  CopeNumber <- unique(samples$CopeNumber)
  type <- unique(samples$type)
  iter <- unique(samples$iter)
  n_sub <- unique(samples$n_sub)

  .get_glm(
    files = files,
    Task = Task,
    CopeNumber = CopeNumber,
    type = type,
    iter = iter,
    n_sub = n_sub
  )
}

get_glm_pop <- function(test, contrasts = NULL) {
  if ("MSMALL" %in% names(test)) {
    files <- test |>
      tidyr::pivot_longer(
        c(MSMALL, SURFACE, VOL),
        names_to = "type",
        values_to = "src"
      ) |>
      dplyr::select(-tar_group) |>
      na.omit() |>
      dplyr::semi_join(
        contrasts,
        by = dplyr::join_by(Task, CopeNumber, ContrastName)
      )
  } else {
    files <- test |>
      tidyr::pivot_longer(
        UKB,
        names_to = "type",
        values_to = "src"
      ) |>
      na.omit()
  }

  CopeNumber <- unique(files$CopeNumber)
  Task <- unique(files$Task)

  out <- tibble::tibble()
  for (type in unique(files$type)) {
    files_by_type <- dplyr::filter(files, type == .env$type)

    out <- dplyr::bind_rows(
      out,
      .get_glm(
        files = files_by_type,
        Task = Task,
        CopeNumber = CopeNumber,
        type = type,
        iter = 0,
        n_sub = nrow(files_by_type)
      )
    )
  }
  out
}

get_glm_pop_pre_ukb <- function(
  src = "data-raw/ukb_mni_all_mean.nii.gz",
  src_sd = "data-raw/ukb_mni_all_sd.nii.gz"
) {
  n_sub <- 38347
  iter <- 0
  Task <- "EMOTION"
  CopeNumber <- 5
  type <- "UKB"

  avgs <- src |>
    to_tbl() |>
    mask() |>
    mask_gray() |>
    mask_atlas() |>
    magrittr::use_series("value")

  sds <- src_sd |>
    to_tbl() |>
    mask() |>
    mask_gray() |>
    mask_atlas() |>
    magrittr::use_series("value")
  out_file <- fs::path(
    Sys.getenv("NIIDIR"),
    glue::glue(
      "nsub-{n_sub}_iter-{iter}_task-{Task}_cope-{CopeNumber}_type-{type}_glm.parquet"
    )
  )
  tibble::tibble(pe = avgs, sigma = sds) |>
    dplyr::mutate(
      index = 1:dplyr::n(),
      z = (sqrt(n_sub) * avgs / sds) |> # t_stat -> z
        pt(n_sub - 1, log.p = TRUE, lower.tail = FALSE) |>
        qnorm(log.p = TRUE, lower.tail = FALSE)
    ) |>
    arrow::write_parquet(out_file, compression = "zstd")

  tibble::tibble(
    Task = Task,
    CopeNumber = CopeNumber,
    type = type,
    iter = iter,
    n_sub = n_sub,
    glm = out_file
  )
}

.cor_one_study <- function(glm.study, glm.gold, method = "spearman") {
  dplyr::bind_rows(
    list(
      study = arrow::read_parquet(glm.study),
      gold = arrow::read_parquet(glm.gold)
    ),
    .id = "sim"
  ) |>
    dplyr::mutate(cohens_d = pe / sigma) |>
    dplyr::select(-tidyselect::any_of(c("z", "pe", "sigma", "n_sub"))) |>
    tidyr::pivot_wider(names_from = sim, values_from = cohens_d) |>
    dplyr::filter(!is.na(study), !is.na(gold)) |>
    dplyr::summarise(
      rho = cor(study, gold, method = method)
    ) |>
    magrittr::use_series(rho)
}


make_data_topo_gold_to_study2 <- function(glm, glm_pop, method = "spearman") {
  dplyr::left_join(
    glm,
    dplyr::select(glm_pop, -iter, -n_sub),
    by = dplyr::join_by(Task, CopeNumber, type),
    suffix = c(".study", ".gold")
  ) |>
    dplyr::mutate(
      rho = purrr::map2_dbl(
        glm.study,
        glm.gold,
        .cor_one_study,
        method = method
      )
    ) |>
    dplyr::select(-tidyselect::starts_with("glm"))
}

# make_data_topo_gold_to_study <- function(tfce, tfce_pop, method = "spearman") {
#   tfce_pop <- tfce_pop |>
#     dplyr::semi_join(tfce, by = c("ContrastName"))
#   checkmate::assert_data_frame(tfce, nrows = 1)
#   checkmate::assert_data_frame(tfce_pop, nrows = 1)
#
#
#   study <- get_pairs(fs::path(storage_dir, fs::path_file(tfce$tstat)), tfce$n_sub) |>
#     mask() |>
#     mask_gray() |>
#     dplyr::mutate(study = cope / sigma * correct_d(tfce$n_sub)) |>
#     dplyr::select(x, y, z, study)
#
#   test <- get_pairs(fs::path(storage_dir, fs::path_file(tfce_pop$tstat)), tfce_pop$n_sub) |>
#     mask() |>
#     mask_gray() |>
#     dplyr::mutate(test = cope / sigma * correct_d(tfce_pop$n_sub)) |>
#     dplyr::select(x, y, z, test)
#
#   tfce |>
#     dplyr::select(Task, CopeNumber, ContrastName, iter, n_sub) |>
#     dplyr::bind_cols(
#       dplyr::left_join(study, test, by = c("x", "y", "z")) |>
#         dplyr::summarise(rho = cor(study, test, method = .env$method, use = "complete.obs"))
#     ) |>
#     dplyr::mutate(method = .env$method)
# }

make_data_topo_study_to_study <- function(glm, times = 2000) {
  glm |>
    dplyr::mutate(
      glm = stringr::str_remove(glm, "/dcl01/smart/data/psadil/meta/"),
      data = purrr::map(
        glm,
        ~ arrow::open_dataset(.x) |>
          dplyr::mutate(d = pe / sigma) |>
          dplyr::select(d, index) |>
          dplyr::collect()
      )
    ) |>
    dplyr::select(-glm) |>
    tidyr::unnest(data) |>
    tidyr::pivot_wider(names_from = iter, values_from = d) |>
    dplyr::group_nest(Task, CopeNumber, type, n_sub, tar_group) |>
    tidyr::crossing(method = c("consistency", "agreement")) |>
    dplyr::mutate(
      fit = purrr::map2(
        data,
        method,
        ~ .x |>
          dplyr::select(-index) |>
          irr::icc(model = "two", type = .y)
      ),
      estimate = purrr::map_dbl(fit, ~ .x$value),
      lower = purrr::map_dbl(fit, ~ .x$lbound),
      upper = purrr::map_dbl(fit, ~ .x$ubound)
    ) |>
    dplyr::select(-fit, -data)
}


.cor_one_study_bynetwork <- function(
  glm.study,
  glm.gold,
  at
) {
  dplyr::bind_rows(
    list(
      study = arrow::read_parquet(glm.study),
      gold = arrow::read_parquet(glm.gold)
    ),
    .id = "sim"
  ) |>
    dplyr::mutate(cohens_d = pe / sigma) |>
    dplyr::select(-tidyselect::any_of(c("z", "pe", "sigma", "n_sub"))) |>
    tidyr::pivot_wider(names_from = sim, values_from = cohens_d) |>
    dplyr::filter(!is.na(study), !is.na(gold)) |>
    dplyr::left_join(at, by = dplyr::join_by(index)) |>
    dplyr::summarise(
      rho = cor(study, gold, method = "spearman"),
      .by = `Network Name`
    )
}

make_data_topo_gold_to_study_bynetwork <- function(glm2, glm_pop2, at) {
  mapping <- to_tbl(MNITemplate::getMNIPath("Brain", res = "2mm")) |>
    mask() |>
    mask_gray() |>
    mask_atlas() |>
    dplyr::mutate(index = 1:dplyr::n()) |>
    dplyr::select(-value)

  at2 <- mapping |>
    dplyr::left_join(at, by = dplyr::join_by(x, y, z)) |>
    dplyr::mutate(
      `Network Name` = dplyr::if_else(
        is.na(`Network Name`) & !is.na(label),
        "subcortical",
        `Network Name`
      )
    ) |>
    dplyr::select(`Network Name`, index) |>
    na.omit()

  dplyr::left_join(
    glm2,
    dplyr::select(glm_pop2, -iter, -n_sub),
    by = dplyr::join_by(Task, CopeNumber, type),
    suffix = c(".study", ".gold")
  ) |>
    dplyr::filter(type == "VOL") |>
    dplyr::mutate(dplyr::across(
      tidyselect::starts_with("glm."),
      ~ stringr::str_remove(.x, "/dcl01/smart/data/psadil/meta/")
    )) |>
    dplyr::mutate(
      rho = purrr::map2(
        glm.study,
        glm.gold,
        ~ .cor_one_study_bynetwork(.x, .y, at = at2),
      )
    ) |>
    dplyr::select(-tidyselect::starts_with("glm")) |>
    tidyr::unnest(rho) |>
    dplyr::mutate(
      Task = stringr::str_to_lower(Task),
      `N Sub` = glue::glue("N Sub: {n_sub}"),
      `N Sub` = factor(
        `N Sub`,
        levels = c(
          "N Sub: 20",
          "N Sub: 40",
          "N Sub: 60",
          "N Sub: 80",
          "N Sub: 100"
        )
      )
    )
}
