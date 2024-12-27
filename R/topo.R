
.get_d <- function(
    files,
    Task,
    CopeNumber,
    type,
    iter,
    n_sub
){
  pe_stem <- fs::path(
    Sys.getenv("NIIDIR"), 
    glue::glue("nsub-{n_sub}_iter-{iter}_task-{Task}_cope-{CopeNumber}_pe")
  )
  sigma_stem <- fs::path(
    Sys.getenv("NIIDIR"),
    glue::glue("nsub-{n_sub}_iter-{iter}_task-{Task}_cope-{CopeNumber}_sigma")
  )
  z_stem <- fs::path(
    Sys.getenv("NIIDIR"),
    glue::glue("nsub-{n_sub}_iter-{iter}_task-{Task}_cope-{CopeNumber}_z")
  )
  
  if (type %in% c("MSMALL", "SURFACE")){
    xiis <- purrr::map(files$src, ciftiTools::read_cifti) |>
      ciftiTools::merge_xifti(xifti_list = _) 
    avgs <- ciftiTools::apply_xifti(xiis, margin=1, FUN=mean) 
    sds <- ciftiTools::apply_xifti(xiis, margin=1, FUN=sd) 
    t_stat <- sqrt(n_sub) * avgs / sds
    z_stat <- ciftiTools::transform_xifti(
      t_stat, 
      FUN=\(x) pt(x, n_sub-1, log.p = TRUE, lower.tail = TRUE)
    ) |>
      ciftiTools::transform_xifti(
        FUN=\(x) qnorm(x, log.p = TRUE, lower.tail = FALSE)
      )
    
    extension <- ".dscalar.nii"
    ciftiTools::write_cifti(avgs, glue::glue("{pe_stem}{extension}"))
    ciftiTools::write_cifti(sds, glue::glue("{sigma_stem}{extension}"))
    ciftiTools::write_cifti(z_stat, glue::glue("{z_stem}{extension}"))
    
  }else{
    copes <- RNifti::readNifti(files$src) |> simplify2array()
    stopifnot(length(dim(copes)) == 4)
    avgs_arr <- apply(copes, 1:3, mean)
    sds_arr <- apply(copes, 1:3, sd)
    t_stat <- sqrt(n_sub) * avgs_arr / sds_arr
    z_stat_arr <- apply(
      t_stat, 1:3,
      FUN = function(x) {
        pt(x, n - 1, log.p = TRUE, lower.tail = FALSE) |>
          qnorm(log.p = TRUE, lower.tail = FALSE)
      }
    )
    z_stat_arr[is.na(z_stat_arr)] <- 0
    
    z_stat <- neurobase::niftiarr(neurobase::readnii(cope_files[[1]]), z_stat_arr)
    avgs <- neurobase::niftiarr(neurobase::readnii(cope_files[[1]]), avgs_arr)
    sds <- neurobase::niftiarr(neurobase::readnii(cope_files[[1]]), sds_arr)
    
    extension <- ".nii.gz"
    RNifti::writeNifti(avgs, glue::glue("{pe_stem}{extension}"))
    RNifti::writeNifti(sds, glue::glue("{sigma_stem}{extension}"))
    RNifti::writeNifti(z_stat, glue::glue("{z_stem}{extension}"))
  }
  
  tibble::tibble(
    Task = Task,
    CopeNumber = CopeNumber,
    type = type,
    iter = iter,
    n_sub = n_sub,
    pe = glue::glue("{pe_stem}{extension}"),
    sigma = glue::glue("{sigma_stem}{extension}"),
    zstat = glue::glue("{z_stem}{extension}")
  )
}


get_d <- function(test, hcp_samples){
  
  files <- test |> 
    tidyr::pivot_longer(
      c(MSMALL, SURFACE, VOL), 
      names_to = "type", 
      values_to = "src") |> 
    dplyr::select(-tar_group) |>
    na.omit() |>
    dplyr::right_join(
      hcp_samples, 
      by = dplyr::join_by(Task, sub, CopeNumber, type))
  
  Task <- unique(hcp_samples$Task)
  CopeNumber <- unique(hcp_samples$CopeNumber)
  type <- unique(hcp_samples$type)
  iter <- unique(hcp_samples$iter)
  n_sub <- unique(hcp_samples$n_sub)
  
  .get_d(
    files=files, 
    Task=Task, 
    CopeNumber = CopeNumber, 
    type=type, 
    iter=iter, 
    n_sub=n_sub)
}

get_d_pop <- function(test, contrasts){
  
  files <- test |> 
    tidyr::pivot_longer(
      c(MSMALL, SURFACE, VOL), 
      names_to = "type", 
      values_to = "src") |> 
    dplyr::select(-tar_group) |>
    na.omit() |>
    dplyr::semi_join(
      contrasts, 
      by = dplyr::join_by(Task, CopeNumber, ContrastName)
    )
  
  CopeNumber <- unique(contrasts$CopeNumber)
  Task <- unique(contrasts$Task)
  iter <- 0
  
  out <- tibble::tibble()
  for (type in unique(files$type)){
    files_by_type <- dplyr::filter(files, type==.env$type)
    
    out <- dplyr::bind_rows(
      out,
      .get_d(
        files=files_by_type, 
        Task=Task, 
        CopeNumber=CopeNumber, 
        type=type, 
        iter=iter, 
        n_sub=nrow(files_by_type))
    )
  }
  out
}



