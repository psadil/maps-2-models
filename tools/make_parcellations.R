# File locations
dst <- "data-raw/parc.rda"

library(ciftiTools)

# from https://github.com/mandymejia/ciftiTools/blob/a1dfa46dee850111bd6938859980da458977c26f/data-raw/make_data_and_sysdata.R

n_parcels <- c(200, 400, 600, 800, 1000)


# Compress the parcellations
parc_name <- glue::glue("data-raw/Parcellations/HCP/fslr32k/cifti/Schaefer2018_{n_parcels}Parcels_7Networks_order.dlabel.nii")
parc <- lapply(parc_name, read_xifti)
names(parc) <- gsub(".dlabel.nii", "", parc_name, fixed=TRUE) |>
  stringr::str_remove("data-raw/Parcellations/HCP/fslr32k/cifti/")
parc <- lapply(parc, function(y){list(
  map = as.matrix(y),
  col = y$meta$cifti$labels$parcels[c("Red", "Green", "Blue")]
)})


save(parc, file=dst, compress='xz')
