get_stats_dirs <- function(root){
  fs::dir_ls(root, recurse = TRUE, glob = "*stats*", type = "directory") %>%
    fs::path_dir()
}

prep_cope_tbl <- function(cope_files, n, study, iter=1){
  copes <- sample(cope_files, n)
  tibble::tibble(copes = copes, n_sub=n, study=study, iter=iter)
}


calc_z <- function(cope_files){
  copes <- purrr::map(cope_files, ~neurobase::readnii(.x)@.Data) %>%
    simplify2array()
  stopifnot(length(dim(copes)) == 4)
  n <- dim(copes)[4]
  means <- apply(copes, 1:3, mean)
  sds <- apply(copes, 1:3, sd)
  t_stat <- sqrt(n) * means / sds
  z_stat <- apply(t_stat, 1:3, FUN = function(x) qnorm(pt(x, n-1, log.p = TRUE), log.p = TRUE))
  z_stat[is.na(z_stat)] <- 0
  
  neurobase::niftiarr(neurobase::readnii(cope_files[[1]]), z_stat)
}


calc_clusters <- function(cope_files, pthresh = 0.05){
  # TODO: standardize to mni
  
  z_stat <- calc_z(cope_files$copes)
  z_file <- neurobase::writenii(z_stat, fs::file_temp())
  
  #' /vols/Data/ukbiobank/local/bb_FSL/bin/cluster 
  #' -i thresh_zstat1 
  #' -c stats/cope1 
  #' -t 2.3 
  #' -p 0.05 
  #' -d 0.388263 
  #' -x reg/example_func2highres.mat 
  #' --warpvol=reg/highres2standard_warp 
  #' --stdvol=reg/standard 
  #' --mm
  #' --volume=119820 
  #' --othresh=thresh_zstat1 
  #' -o cluster_mask_zsqtat1 
  #' --connectivity=26  
  #' --olmax=lmax_zstat1.txt 
  #' --scalarname=Z > cluster_zstat1.txt
  
  clusters <- fslr::fsl_cluster(
    file = z_file,
    threshold = 2.3,
    pthresh = pthresh,
    smooth_est = 0.388263,
    opts = glue::glue("--volume={119820}"),
    mm = TRUE,
    connectivity = 26
  )
  
  fslr::read_cluster_table(clusters$olmax) %>%
    tibble::as_tibble() %>%
    dplyr::rename(index = `Cluster Index`) %>%
    dplyr::mutate(
      n_sub = unique(cope_files$n_sub),
      study = unique(cope_files$study),
      iter = unique(cope_files$iter))
}


write_study <- function(cl, study, n_sub, iter){
  file <- fs::file_temp()
  
  index <- if (dplyr::n_distinct(cl$index) > 1) "all" else unique(cl$index)
  
  header <- c(
    glue::glue("// index {index}, study {study}, iter {iter}"), 
    glue::glue("// Subjects={n_sub}"))
  
  readr::write_lines(header, file)
  cl %>%
    dplyr::select(x, y, z) %>%
    readr::write_tsv(file, append = TRUE) 
  
  # add empty line, for separating studies
  readr::write_lines("", file, append = TRUE)
  
  return(file)
}

write_all_studies <- function(study_files, out_file){
  readr::write_lines("// Reference=MNI", out_file)
  readr::write_lines(readr::read_lines(study_files), out_file, append = TRUE)
  invisible(out_file)
}

do_ale <- function(
  clusters, 
  perm=1000, 
  p=0.01, 
  clust=0.01){
  
  cl <- clusters %>%
    dplyr::group_by(study, n_sub, iter) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      cl_files = purrr::pmap_chr(
        list(cl=data, study=study, n_sub=n_sub, iter=iter), 
        write_study))
  
  iter <- unique(clusters$iter)
  n_sub <- unique(clusters$n_sub)
  index <- if (dplyr::n_distinct(clusters$index) > 1) "all" else unique(clusters$index)
  
  foci_file <- fs::path_temp(glue::glue("foci_nsub-{n_sub}_iter-{iter}_index-{index}")) %>%
    fs::path_abs()
  
  write_all_studies(cl$cl_files, foci_file)
  
  system2(
    "java",
    # args = glue::glue("-cp data-raw/GingerALE.jar org.brainmap.meta.getALE2 {foci_file} -mask=/home/psadil/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz -p={p} -perm={perm} -clust={clust}")
    # args = glue::glue("-cp data-raw/GingerALE.jar org.brainmap.meta.getALE2 {foci_file} -p={p} -mask=MNI152_wb -minVol=20")
    args = glue::glue("-cp data-raw/GingerALE.jar org.brainmap.meta.getALE2 {foci_file} -p={p} -mask=MNI152_wb -clust={clust} -perm={perm}")
  )
  
  fname <- fs::path_file(foci_file)
  foci_dir <- fs::path_dir(foci_file)
  niis <- fs::dir_ls(foci_dir, glob = glue::glue("*{fname}*nii"))
  
  # out_dir <- (foci_stem %||% foci_file) %>%
  #   fs::path_dir() %>%
  #   fs::path_abs()
  out_dir <- here::here("data-raw", "niis")
  
  purrr::walk(niis, R.utils::gzip, overwrite=TRUE)
  # gzs <- fs::dir_ls(foci_dir, glob = glue::glue("*{fname}*nii.gz"))
  
  # out_stem <- foci_file %>%
  #   fs::path_file() %>%
  #   fs::path_ext_remove()
  
  # out_files <- stringr::str_replace(gzs, fname, out_stem)
  # out <- fs::file_move(gzs, out_files) %>%
  #   fs::file_move(out_dir)
  
  if(!fs::dir_exists(out_dir)) fs::dir_create(out_dir)
  
  out <- fs::dir_ls(foci_dir, glob = glue::glue("*{fname}*nii.gz")) %>%
    fs::file_move(out_dir)
  
  out
}

apply_reg_cope <- function(
  feat_dir, 
  outfile = fs::file_temp(),
  reffile = fs::path(Sys.getenv("FSLDIR"), "data", "standard","MNI152_T1_2mm_brain.nii.gz")){
  
  # Option -s ( standard image ) selected with  argument "/vols/Data/ukbiobank/local/bb_FSL/data/standard/MNI152_T1_2mm_brain"
  # convert_xfm -omat example_func2standard.mat -concat highres2standard.mat example_func2highres.mat
  # convertwarp --ref=standard --premat=example_func2highres.mat --warp1=highres2standard_warp --out=example_func2standard_warp
  # applywarp --ref=standard --in=example_func --out=example_func2standard --warp=example_func2standard_warp
  
  warp <- fs::path(feat_dir, "reg", "example_func2standard_warp")
  cope5 <- fs::path(feat_dir, "stats", "cope5")
  
  outfile <- fs::path_abs(outfile)
  
  fslr::fsl_applywarp(
    infile = cope5,
    outfile = outfile,
    reffile = reffile,
    warpfile = warp,
    retimg = FALSE)
  
  invisible(fs::path_ext_set(outfile, ".nii.gz"))
}

