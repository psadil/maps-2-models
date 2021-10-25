library(tidyverse)

tar_load(z_img)

d <- cope_files2 |>
  select(-tar_group) |>
  filter(n_study == 5 ) |>
  rename(z = copes) 

d <- z_img |>
  filter(iter == 1) |>
  mutate(exp = glue::glue("study-{study}_nsub-{n_sub}_nstudy-{n_study}_iter-{iter}")) 

reticulate::use_condaenv("meta")
reticulate::source_python("python/ale.py")


mask_file <- "/home/psadil/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz"

x <- do_ale(d[1:2,], here::here(), here::here(), "tmp", mask_file)
z <- neurobase::niftiarr(neurobase::readnii(mask_file), neurobase::readnii(x))


out <- neurobase::img_indices(z, mask = mask, add_values = TRUE) |>
  tibble::as_tibble()



