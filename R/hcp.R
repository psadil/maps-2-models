not_avail <- function() {
  c(
    110613, 113417, 113821, 120010, 121719, 130518, 139637, 143830, 146836,
    168139, 175035, 176239, 185038, 189652, 199958, 201515, 202820, 385046,
    401422, 415837, 433839, 462139, 465852, 469961, 644246, 656657, 688569,
    723141, 767464, 872764, 943862, 965367, 969476, 987983, 994273, 433839,
    # https://wiki.humanconnectome.org/docs/HCP%20Data%20Release%20Updates%20Known%20Issues%20and%20Planned%20fixes.html
    # 3T Functional Preprocessing Error of all 3T “RL” fMRI runs in 25 Subjects
    103010, 113417, 116423, 120010, 121719, 127226, 130114, 143830, 169040,
    185038, 189652, 202820, 204218, 329844, 385046, 401422, 462139, 469961,
    644246, 688569, 723141, 908860, 943862, 969476, 971160, 
    196952, 748662, 809252, 144428, 186545, 192237, 223929, 320826, 644044,
    822244, 870861, 947668, 
    # Subjects without Field Maps for Structural scans
    102614, 111009, 111514, 115017, 121416, 130821, 138332, 179952, 299760,
    300618, 392750, 406432, 429040, 633847, 662551, 679770, 688569, 693461,
    815247, 
    # duplicated name
    142626
  )
}

cor_w_pop_by_region <- function(tfce, tfce_pop, at, storage_dir, method = "spearman") {
  tfce <- tfce |>
    dplyr::filter(stringr::str_detect(tfce_corrp_tstat, glue::glue("flags-{ContrastName}_tfce")))
  tfce_pop <- tfce_pop |>
    dplyr::semi_join(tfce, by = c("ContrastName"))
  checkmate::assert_data_frame(tfce, nrows = 1)
  checkmate::assert_data_frame(tfce_pop, nrows = 1)
  
  study <- get_pairs(fs::path(storage_dir, fs::path_file(tfce$tstat)), tfce$n_sub) |>
    mask() |>
    mask_gray() |>
    dplyr::mutate(study = cope / sigma * correct_d(tfce$n_sub)) |>
    dplyr::select(x, y, z, study)
  
  test <- get_pairs(fs::path(storage_dir, fs::path_file(tfce_pop$tstat)), tfce_pop$n_sub) |>
    mask() |>
    mask_gray() |>
    dplyr::mutate(test = cope / sigma * correct_d(tfce$n_sub)) |>
    dplyr::select(x, y, z, test)
  
  dplyr::left_join(study, test, by = c("x", "y", "z")) |>
    dplyr::left_join(at, by = c("x", "y", "z")) |>
    dplyr::filter(!is.na(label)) |>
    dplyr::group_by(label, `Label Name`, `Network Name`, `Full component name`, hemi) |>
    dplyr::summarise(rho = cor(study, test, method = .env$method), .groups = "drop") |>
    dplyr::bind_cols(
      tfce |>
        dplyr::select(Task, CopeNumber, ContrastName, iter, n_sub)
    ) |>
    dplyr::mutate(method = .env$method)
}


get_pairs <- function(tstat, n_sub) {
  cope <- to_tbl(stringr::str_replace(tstat, "_tstat1", "_glm_cope"), measure = "cope")
  varcope <- to_tbl(stringr::str_replace(tstat, "_tstat1", "_glm_varcope"), measure = "varcope") |>
    dplyr::mutate(sigma = sqrt(varcope) * sqrt(n_sub)) |>
    dplyr::select(-varcope)
  
  dplyr::left_join(cope, varcope, by = c("x", "y", "z")) |>
    mask() 
}


get_pop_d <- function(tfce_pop) {
  tfce_pop |>
    dplyr::select(tstat, Task, CopeNumber, ContrastName, n_sub) |>
    dplyr::mutate(data = purrr::map2(tstat, n_sub, get_pairs)) |>
    dplyr::select(-tstat) |>
    tidyr::unnest(data) |>
    mask() |>
    dplyr::mutate(
      d = cut(
        abs(cope / sigma * correct_d(n_sub)),
        right = FALSE,
        ordered = TRUE,
        breaks = c(0, .2, .5, .8, Inf),
        labels = c("0", "small", "medium", "large")
      )
    )
}

get_hcp_copes <- function(
    contrasts, 
    msmall_feat_file="data-raw/msmall_feat",
    matching_only = TRUE){
  
  
  hcp_copes <- tibble::tibble(msmall_feat = fs::path(readr::read_lines(msmall_feat_file))) |>
    dplyr::mutate(
      Task = stringr::str_extract(
        msmall_feat, 
        "WM|GAMBLING|MOTOR|LANGUAGE|SOCIAL|RELATIONAL|EMOTION"),
      sub = stringr::str_extract(msmall_feat, "[[:digit:]]{6}")) |>
    dplyr::semi_join(get_hcp_twins(), by = dplyr::join_by(sub)) |>
    dplyr::inner_join(contrasts, by = dplyr::join_by(Task)) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      MSMALL = fs::path_join(c(msmall_feat, "GrayordinatesStats", glue::glue("cope{CopeNumber}.feat"))),
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-msmall_feat) |>
    dplyr::mutate(
      VOL = stringr::str_replace(MSMALL, "_MSMAll.feat/GrayordinatesStats/", "vol.feat/"),
      VOL = stringr::str_c(VOL, "/stats"),
      SURFACE = stringr::str_replace(MSMALL, "_MSMAll.feat", ".feat")
    ) |>
    tidyr::pivot_longer(c(MSMALL, VOL, SURFACE)) |>
    dplyr::mutate(
      copefile = dplyr::case_match(
        name,
        c("MSMALL", "SURFACE") ~ "cope1.dtseries.nii",
        "VOL" ~ "cope1.nii.gz"
      ),
      value = stringr::str_c(value, copefile, sep="/"),
      exists = fs::file_exists(value)) |>
    dplyr::filter(exists) |>
    dplyr::select(-exists, -copefile) |>
    tidyr::pivot_wider() 
  
  if (matching_only) {
    hcp_copes |>
      na.omit()
  }else{
    hcp_copes 
  }
  
}

get_hcp_twins <- function(unrestricted="data-raw/hcp/restricted.csv"){
  readr::read_csv(
    unrestricted,
    col_select = c(Subject, ZygositySR, ZygosityGT, Mother_ID),
    show_col_types = FALSE) |>
    dplyr::rename(sub=Subject) |>
    dplyr::mutate(
      ZygositySR = dplyr::if_else(ZygositySR=="NotMZ", "DZ", ZygositySR),
      zygosity = dplyr::if_else(is.na(ZygosityGT), ZygositySR, ZygosityGT),
      twin_group = dplyr::case_when(
        zygosity == "NotTwin" ~ sub,
        stringr::str_detect(zygosity, "Z") ~ Mother_ID,
        TRUE ~ NA
      )
    ) |>
    dplyr::select(sub, twin_group) |>
    dplyr::filter(!sub %in% not_avail()) |>
    dplyr::arrange(twin_group) |>
    dplyr::group_by(twin_group) |>
    dplyr::slice_head(n=1) |>
    dplyr::ungroup() |>
    dplyr::mutate(sub=as.character(sub)) 
}

get_hcp_contrasts <- function(src=here::here("data-raw", "hcp", "contrasts.csv")){
  readr::read_csv(src, show_col_types = FALSE) |>
    dplyr::filter(
      (Task == "EMOTION" & CopeNumber %in% c(3)) |
        (Task == "WM" & CopeNumber %in% c(11)) |
        (Task == "GAMBLING" & CopeNumber %in% c(6)) |
        (Task == "MOTOR" & CopeNumber %in% c(7)) |
        (Task == "LANGUAGE" & CopeNumber %in% c(3)) |
        (Task == "SOCIAL" & CopeNumber %in% c(3)) |
        (Task == "RELATIONAL" & CopeNumber %in% c(3))
    )
}

sample_hcp <- function(
    test, 
    n_iter, 
    n_subs, 
    types = c("VOL","MSMALL","SURFACE")){
  
  test |>
    dplyr::select(-tar_group, -ContrastName) |>
    tidyr::pivot_longer(c(SURFACE, MSMALL, VOL), names_to = "type") |>
    dplyr::filter(type %in% c(types)) |>
    na.omit() |>
    tidyr::crossing(
      iter = seq_len(n_iter),
      n_sub = n_subs) |>
    dplyr::group_nest(Task, CopeNumber, type, iter, n_sub) |>
    dplyr::mutate(
      data = purrr::map2(data, n_sub, ~dplyr::sample_n(.x, size=.y, replace=TRUE))
    ) |>
    tidyr::unnest(data) |>
    dplyr::select(-value)
  
}
