#' TODO: distances in surface

#' make gold-standard
#' do the below for gold 
#' 
#' ??? which "corrected areas" should I used to calculate distances from gold-standard -> study sample ??

#' -cifti-create-dense-from-template
dr <- "/Users/psadil/git/manuscripts/maps-to-models/meta/data-raw"
stem0 <- "nsub-20_iter-85_task-WM_cope-11_type-MSMALL"
stem <- glue::glue("{dr}/palm/{stem0}")
fwe <- glue::glue("{dr}/palm/fwep.dscalar.nii")
volume <- glue::glue("{dr}/hcp-niis-ptfce/{stem0}_tfce_tstat_fwep_c1.nii")
lh <- glue::glue("{dr}/hcp-niis-ptfce/{stem0}_L_tfce_tstat_fwep_c1.gii")
rh <- glue::glue("{dr}/hcp-niis-ptfce/{stem0}_R_tfce_tstat_fwep_c1.gii")


dense <- ciftiTools::run_wb_cmd(
  glue::glue("-cifti-create-dense-from-template {stem}.dtseries.nii {fwe} -volume-all {volume} -metric CORTEX_LEFT {lh} -metric CORTEX_RIGHT {rh}"),
  intern = TRUE
)

volume2 <- glue::glue("{dr}/hcp-niis-ptfce/{stem0}_tfce_tstat_c1.nii")
lh2 <- glue::glue("{dr}/hcp-niis-ptfce/{stem0}_L_tfce_tstat_c1.gii")
rh2 <- glue::glue("{dr}/hcp-niis-ptfce/{stem0}_R_tfce_tstat_c1.gii")
tstat <- glue::glue("{dr}/palm/tstat.dscalar.nii")

dense2 <- ciftiTools::run_wb_cmd(
  glue::glue("-cifti-create-dense-from-template {stem}.dtseries.nii {tstat} -volume-all {volume2} -metric CORTEX_LEFT {lh2} -metric CORTEX_RIGHT {rh2}"),
  intern = TRUE
)

#' threshold tfce_tstat_c1 by tfce_tstat_fwep_c1 (neg).
fwe_cii <- ciftiTools::read_cifti(fwe)
mask <- fwe_cii > -log(0.05)
tstat_cii <- ciftiTools::read_cifti(tstat)
masked <- tstat_cii * mask
masked_file <- glue::glue("{dr}/palm/masked.dscalar.nii")
ciftiTools::write_cifti(masked, masked_file)

#' find extrema with -cifti-extrema
extrema <- glue::glue("{dr}/palm/extrema.dscalar.nii")
l_mid <- glue::glue("{stem}_L_midthickness.32k_fs_LR.gii")
r_mid <- glue::glue("{stem}_R_midthickness.32k_fs_LR.gii")

ciftiTools::run_wb_cmd(
  glue::glue("-cifti-extrema {masked_file} 10 10 COLUMN {extrema} -left-surface {l_mid} -right-surface {r_mid}"),
  intern = FALSE
)


#' in R, create table of indices/extrema, then filter and store
ee <- ciftiTools::read_cifti(extrema)
left_max <- which(ee$data$cortex_left==1)

#' -surface-geodesic-distance to make dconn files (one file per extrema)
#' wb_command -surface-geodesic-distance nsub-20_iter-85_task-WM_cope-11_type-MSMALL_L_midthickness.32k_fs_LR.gii 973 l_dist.shape.gii  -limit 20 -corrected-areas nsub-20_iter-85_task-WM_cope-11_type-MSMALL_L_area.32k_fs_LR.gii
#' 
#' load dconns, filter such that to/from are both extrema

