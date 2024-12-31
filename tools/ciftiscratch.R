library(dplyr)

# testing out what happens with surface-average

gii0 <- "/Users/psadil/Desktop/ciftitmp/147737.L.midthickness_MSMAll.32k_fs_LR.surf.gii"
gii1 <- "/Users/psadil/Desktop/ciftitmp/165840.L.midthickness_MSMAll.32k_fs_LR.surf.gii"

surf0 <- freesurferformats::read.fs.surface(gii0)
surf1 <- freesurferformats::read.fs.surface(gii1)

((surf0$vertices + surf1$vertices) / 2) |> 
  head()

merged <- freesurferformats::read.fs.surface("/Users/psadil/Desktop/ciftitmp/merged.surf.gii")

merged$vertices |> head()

merged$faces |> head()
surf0$faces |> head()
surf1$faces |> head()

all.equal(surf0$faces, surf1$faces)

# confirm what happens with extrema
left_surface <- "data-raw/palm/nsub-379_iter-0_task-WM_cope-11_type-MSMALL_L_midthickness.32k_fs_LR.gii"
right_surface <- "data-raw/palm/nsub-379_iter-0_task-WM_cope-11_type-MSMALL_R_midthickness.32k_fs_LR.gii"
gii <- freesurferformats::read.fs.surface(left_surface)

faces <- as_tibble(gii$faces)

connected_to_1 <- faces |>
  filter(V1==1 | V2==1 | V3==1) |>
  tidyr::pivot_longer(everything()) |>
  distinct(value) |>
  magrittr::use_series("value") |>
  sort()

xii <- ciftiTools::read_cifti(cohens_gold) |>
  ciftiTools::move_from_mwall()

xii$data$cortex_left[connected_to_1]

new_file <- tempfile(fileext = ".dscalar.nii")
xii$data$cortex_left[1] <- 10
ciftiTools::write_cifti(xii, new_file)

cifti_out <- tempfile(fileext = ".dscalar.nii")

ciftiTools::run_wb_cmd(
  glue::glue("-cifti-extrema {new_file} 1 1 COLUMN {cifti_out} -left-surface {left_surface} -right-surface {right_surface} -only-maxima")
)

xii_max <- ciftiTools::read_cifti(cifti_out) 
# conclusion: cifti-extrema cares only about connectivity, plus the area part
 

# is the geodesic distance really just a search along a weighted graph?

# simple, without connected areas
distances2 <- tempfile(fileext = ".shape.gii")
ciftiTools::run_wb_cmd(
  glue::glue("-surface-geodesic-distance {left_surface} 0 {distances2}"),
)

gii_distances2 <- freesurferformats::read.fs.morph(distances2)
gii_distances2[connected_to_1]

vertices <- as_tibble(gii$vertices)
dist(rbind(gii$vertices[1,],gii$vertices[13,]))
dist(rbind(gii$vertices[1,],gii$vertices[69,]))
# conclusion: geodesic distance is building a path of euclidean distances between vertices

# what happens when using corrected_areas?
distances <- tempfile(fileext = ".shape.gii")
corrected_areas <- "data-raw/palm/nsub-379_iter-0_task-WM_cope-11_type-MSMALL_L_area.32k_fs_LR.gii"
ciftiTools::run_wb_cmd(
  glue::glue("-surface-geodesic-distance {left_surface} 0 {distances} -corrected-areas {corrected_areas}"),
)

areas <- freesurferformats::read.fs.morph(corrected_areas)
areas[connected_to_1]
gii_distances <- freesurferformats::read.fs.morph(distances)
gii_distances[connected_to_1]

# seems like the use of corrected-areas does something complicated with a combination
# of current areas and correction areas:
# https://github.com/Washington-University/workbench/blob/f32180853d5744dc5441df7539bed4a558eada0f/src/Files/GeodesicHelper.cxx#L47



