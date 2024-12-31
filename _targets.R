# TODO: double-check z-stat calculations (dof/cifti writing/header)
# TODO: make permutation test for ciftis

#' TODO: distances in surface
#' -cifti-create-dense-from-template
#' threshold tfce_tstat_c1 by tfce_tstat_fwep_c1 (neg).
#' find extrema with -cifti-extrema
#' -cifti-export-dense-mapping to produce indices
#' in R, create table of indices/extrema, then filter and store
#' -surface-geodesic-distance to make dconn files (one file per extrema)
#' load dconns, filter such that to/from are both extrema

library(targets)
library(tarchetypes)
library(rlang)

source(here::here("R", "spatial.R"))
source(here::here("R", "tfce.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "hcp.R"))
source(here::here("R", "ptfce.R"))
source(here::here("R", "figures.R"))
source(here::here("R", "pairwise.R"))
source(here::here("R", "model.R"))
source(here::here("R", "manuscript.R"))
source(here::here("R", "roi.R"))
source(here::here("R", "topo.R"))
source(here::here("R", "ukb.R"))
source(here::here("R", "cifti.R"))

Sys.setenv(
  NIIDIR = here::here("data-raw", "hcp-niis-ptfce"),
  PALMDIR = here::here("data-raw", "palm"),
  HCPPARQUET = "/Users/psadil/Library/CloudStorage/OneDrive-JohnsHopkins/data/hcp-to-parquet/data/out",
  PALMBIN = "/fastscratch/myscratch/pssadil/PALM/palm"
) # explicitly avoiding tracking this

controller_small <- crew::crew_controller_local(
  name = "small",
  workers = 1
)

# controller_large <- crew::crew_controller_local(
#   name = "large",
#   workers = 5,
# )

# controller_slurm <- crew.cluster::crew_controller_slurm(
#   name = "slurm",
#   workers = 2,
#   script_lines = "source ~/.bashrc; mamba activate meta; ml fsl",
#   slurm_memory_gigabytes_per_cpu=2,
#   slurm_cpus_per_task=1,
#   slurm_time_minutes=10,
#   slurm_partition="shared",
#   command_submit = '/usr/bin/sbatch --constraint="intel"'
# )


targets::tar_option_set(
  format = "qs",
  storage = "worker",
  packages = c("oro.nifti"),
  controller = crew::crew_controller_group(controller_small),
  resources = tar_resources(
    crew = tar_resources_crew(controller = "small")
  )
)


list(
  tar_group_by(
    contrasts,
    get_hcp_contrasts(),
    Task, CopeNumber,
    format = "parquet"
  ),
  tar_target(test, get_hcp_copes(contrasts), format = "parquet"),
  tar_target(test_ukb, get_ukb_copes("data-raw/ukb_copes"), format = "parquet"),
  tar_target(
    test_all, 
    get_hcp_copes(contrasts, matching_only=FALSE), 
    format = "parquet"),
  tar_target(n_sub, c(20, 40, 60, 80, 100)),
  tar_target(iter, seq_len(100)),
  tar_group_by(
    hcp_samples,
    sample_hcp(test=test, n_iter=100, n_subs=c(20, 40, 60, 80, 100)),
    Task, CopeNumber, type, iter, n_sub,
    format = "parquet"
  ),
  tar_group_by(
    ukb_samples,
    sample_ukb(test=test_ukb, n_iter=100, n_subs=c(20, 40, 60, 80, 100, 1000, 10000)),
    Task, CopeNumber, type, iter, n_sub,
    format = "parquet"
  ),
  tar_group_by(
    hcp_samples_all,
    sample_hcp(
      test=test_all, 
      n_iter=100, 
      n_subs=c(20, 40, 60, 80, 100),
      types = c("MSMALL")),
    Task, CopeNumber, type, iter, n_sub,
    format = "parquet"
  ),
  tar_target(
    roi_avg,
    avg_roi(test, n_parcels),
    pattern = cross(test, n_parcels),
    format = "parquet"
  ),
  tar_target(
    rois,
    test_roi(roi_avg, hcp_samples),
    pattern = map(hcp_samples),
    format = "parquet"
  ),
  tar_target(rois_pop, test_roi_pop(roi_avg), format = "parquet"),
  tar_target(
    roi_avg_ukb,
    avg_roi_ukb(test_ukb, n_parcels),
    pattern = cross(test_ukb, n_parcels),
    format = "parquet"
  ),
  tar_target(roi_avg_ukb2, roi_avg_ukb, format = "parquet"),
  tar_target(
    rois_ukb,
    test_roi(roi_avg_ukb2, ukb_samples),
    pattern = map(ukb_samples),
    format = "parquet"
  ),
  tar_target(rois_pop_ukb, test_roi_pop(roi_avg_ukb), format = "parquet"),
  # strategy for palm: create commands, then run with slurm array
  tar_target(tfce, get_tfce_cmd(hcp_samples, test), format = "parquet"),
  tar_target(tfce_pop, get_tfce_pop(test), format = "parquet"),
  # this is back in R/targets (assumes above have been created)
  tar_target(
    study_peaks, 
    get_study_peaks_cifti(tfce), 
    format="parquet",
    pattern = head(map(tfce), 10)),
  tar_target(
    study_to_gold_distances,
    get_cifti_augmented2(study_peaks, tfce_pop), 
    format="parquet"
  ),
  tar_target(at, make_atlas_full()),
  tar_target(n_parcels, c(400)),
  tar_target(
    at_list,
    make_atlas_full(n_parcels = n_parcels),
    pattern = map(n_parcels)
  ),
  tar_target(
    ptfce, 
    do_ptfce2(hcp_samples=hcp_samples, test=test, enhance=TRUE),
    pattern = map(hcp_samples),
    format = "parquet"
  ),
  tar_target(ptfce_pop, do_ptfce_pop(test), format = "parquet"),
  tar_target(
    gold_peaks,
    get_ptfce_maxes_pop(ptfce_pop),
    pattern = map(ptfce_pop),
    format = "parquet"
  ),
  tar_target(do_fwe_correction, c(TRUE, FALSE)),
  tar_target(
    maxes,
    get_ptfce_maxes(row=ptfce, do_fwe_correction=do_fwe_correction),
    pattern = cross(map(ptfce), do_fwe_correction),
    format = "parquet"
  ),
  tar_target(
    augmented,
    augment_distance2(maxes=maxes, gold_peaks=gold_peaks),
    format = "parquet"
  ),
  tar_target(
    tfce_null,
    do_ptfce(
      fs::dir_ls("data-raw/Fake2B"),
      n_sub = n_sub,
      iter = iter,
      storage_dir = here::here("data-raw/Fake2B-ptfce"),
      flags = "2BK-0BK",
      resample = TRUE
    ),
    pattern = cross(n_sub, iter)
  ),
  tar_target(
    tfce_pop_null,
    do_ptfce(
      fs::dir_ls("data-raw/Fake2B"),
      n_sub = length(fs::dir_ls("data-raw/Fake2B")),
      iter = 0,
      storage_dir = here::here("data-raw/Fake2B-ptfce"),
      flags = "2BK-0BK",
      resample = FALSE,
      enhance = FALSE
    )
  ),
  tar_target(
    active0_null,
    tfce_null |>
      dplyr::mutate(tmp = purrr::map(ptfce, get_active_ptfce)),
    pattern = map(tfce_null)
  ),
  tar_target(
    active_null,
    active0_null |>
      dplyr::select(-copes) |>
      tidyr::unnest(tmp) |>
      dplyr::select(-ptfce),
    format = format_arrow_table()
  ),
  tar_target(
    roi_avg_all,
    avg_roi(test_all, n_parcels),
    pattern = cross(test_all, n_parcels),
    format = "parquet"
  ),
  tar_target(
    rois_all,
    test_roi(roi_avg_all, hcp_samples_all),
    pattern = map(hcp_samples_all),
    format = "parquet"
  ),
  tar_target(rois_pop_all, test_roi_pop(roi_avg_all), format = "parquet"),
  tar_target(ContrastNames, contrasts$ContrastName),
  tar_target(
    pairwise,
    cor_pairwise_ptfce(
      tfce, 
      ContrastNames, 
      n_sub, 
      storage_dir = Sys.getenv("NIIDIR"),
      method = "spearman"),
    cross(n_sub, ContrastNames) # ,
    # resources = tar_resources(
    #   crew = tar_resources_crew(controller = "small")
    # )
  ),
  tarchetypes::tar_group_by(
    subtask,
    get_subs_tasks(Sys.getenv("HCPPARQUET"), contrasts),
    sub,
    format = "parquet"
  ),
  tar_target(
    data_topo_sub_to_sub,
    do_loo_cor(
      dataset=Sys.getenv("HCPPARQUET"), 
      subtask=dplyr::select(subtask, sub, task), 
      method="spearman"),
    map(subtask),
    format = "parquet" 
  ),
  tar_target(
    data_peak_study_to_study,
    make_data_peak_study_to_study(
      at = at,
      space = space,
      gold_peaks = gold_peaks
    ),
    format = "parquet"
  ),
  tar_target(
    data_roi_study_to_gold,
    make_data_roi_study_to_gold(
      gold_tested = rois_pop,
      rois_tested = rois,
      at_list=at_list,
      data_topo_gold=data_topo_gold
    ),
    format = "parquet"
  ),
  tar_target(
    data_roi_study_to_study,
    make_data_roi_study_to_study(rois_tested = rois),
    format = "parquet"
  ),
  tar_target(
    data_roi_sub_to_sub,
    make_data_roi_sub_to_sub(rois_pop = rois_pop2),
    format = "parquet"
  ),
  tar_target(
    data_peak_study_to_gold,
    make_data_peak_study_to_gold(
      at = at,
      gold_peaks = gold_peaks,
      space = space,
      data_topo_gold=data_topo_gold
    ),
    format = "parquet"
  ),
  tar_target(
    zstats,
    test |>
      tidyr::unnest(avail) |>
      dplyr::mutate(
        zstat = stringr::str_replace(avail, "/cope1.nii", "/zstat1.nii")
      ) |>
      dplyr::select(Task, CopeNumber, ContrastName, zstat),
    format = "parquet"
  ),
  tar_target(
    sub_peaks,
    dplyr::mutate(
      zstats,
      m = purrr::map(
        .data$zstat,
        ~ get_ptfce_maxes_sub(zstat = .x, cluster_thresh = 0.001, minextent = 0)
      )
    ),
    pattern = map(zstats)
  ),
  tar_target(
    space_sub,
    dplyr::left_join(
      sub_peaks,
      gold_peaks,
      by = c("Task", "CopeNumber", "ContrastName"),
      suffix = c(".study", ".ref")
    ) |>
      dplyr::mutate(
        augmented = purrr::map2(
          m.study, m.ref,
          ~ augment_distance(
            study = .x,
            reference = .y,
            vox_mm = 2
          )
        )
      ) |>
      dplyr::select(augmented, Task, CopeNumber) |>
      tidyr::unnest(augmented)
  ),
  tar_target(
    data_peak_sub_to_sub,
    make_data_peak_sub_to_sub(
      at = at,
      space_sub = space_sub,
      gold_peaks = gold_peaks
    ),
    format = "parquet"
  ),
  tar_target(data_topo_gold, get_pop_d(tfce_pop_fsl)),
  tar_target(
    tfce_fsl,
    tfce |>
      dplyr::mutate(
        tmp = purrr::pmap(
          list(copes=copes, flags=ContrastName, iter=iter, n_sub=n_sub),
          do_tfce2,
          n=1,
          storage_dir = Sys.getenv("NIIDIR"),
        )
      ) |>
      tidyr::unnest(tmp),
    pattern = map(tfce)
  ),
  tar_target(
    data_topo_gold_to_study,
    make_data_topo_gold_to_study(
      tfce_fsl, 
      tfce_pop_fsl, 
      storage_dir = Sys.getenv("NIIDIR"),
      method = "spearman"),
    cross(map(tfce_fsl))
  ),
  tar_target(
    data_topo_study_to_study,
    make_data_topo_study_to_study(pairwise = pairwise),
    format = "parquet"
  ),
  tar_target(
    data_model_gold_gold_to_study,
    make_data_model_gold_gold_to_study(
      dataset_gold = here::here("data-raw/out-perm-gold-cpm-sametest"),
      dataset = here::here("data-raw/out-perm-cpm-sametest")
    ),
    format = "parquet"
  ),
  tar_target(
    data_model_study_to_study,
    make_data_model_study_to_study(
      dataset = here::here("data-raw", "out-perm-cpm-preds-sametest")
    ),
    format = "parquet"
  ),
  tar_target(
    data_model_sub_to_sub,
    make_data_model_sub_to_sub(
      features = here::here("data-raw", "cpm-difumo.parquet")
    ),
    format = "parquet"
  ),
  tar_target(
    pop_cor_region,
    cor_w_pop_by_region(
      tfce = tfce_fsl,
      tfce_pop = tfce_pop_fsl,
      at = at,
      storage_dir = Sys.getenv("NIIDIR"),
      method = "spearman"
    ),
    map(tfce_fsl),
    format = "parquet"
  ),
  tar_target(
    roi,
    make_roi(data_roi_study_to_gold, data_roi_study_to_study, data_roi_sub_to_sub),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    fig_roi,
    make_tikz(p = roi, file = "analyses/figures/roi.tex", width = 3, height = 6),
    format = "file"
  ),
  tar_target(
    prop_active_most_active_roi_ptfce_null,
    make_prop_active_most_active_roi_ptfce_null(
      at_list = at_list,
      active_null = active_null,
      iter = iter,
      gold_tested = gold_tested
    ),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    fig_prop_active_most_active_roi_ptfce_null,
    make_tikz(
      p = prop_active_most_active_roi_ptfce_null,
      file = "analyses/figures/prop-active-most-active-roi-ptfce-null.tex",
      width = 5,
      height = 3
    ),
    format = "file"
  ),
  tar_target(
    prop_active_most_active_roi_ptfce,
    make_prop_active_most_active_roi_ptfce(
      data_roi_study_to_gold = data_roi_study_to_gold
    ),
    packages = c("ggplot2")
  ),
  tar_target(
    fig_prop_active_most_active_roi_ptfce,
    make_tikz(
      p = prop_active_most_active_roi_ptfce,
      file = "analyses/figures/prop-active-most-active-roi-ptfce.tex",
      width = 6,
      height = 6.5
    ),
    format = "file"
  ),
  tar_target(
    peaks,
    make_peaks(
      data_peak_study_to_gold = data_peak_study_to_gold,
      data_peak_study_to_study = data_peak_study_to_study,
      data_peak_sub_to_sub = data_peak_sub_to_sub
    ),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    fig_peaks,
    make_tikz(p = peaks, file = "analyses/figures/peaks.tex", width = 7, height = 3.5),
    format = "file"
  ),
  tar_target(
    peak_bysize,
    make_peak_bysize(space = space, data_topo_gold = data_topo_gold),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    fig_peak_bysize,
    make_tikz(
      p = peak_bysize,
      file = "analyses/figures/peak-bysize.tex",
      width = 5,
      height = 4
    ),
    format = "file"
  ),
  tar_target(
    peak_bynetwork,
    make_peak_bynetwork(space = space, data_topo_gold = data_topo_gold),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    fig_peak_bynetwork,
    make_tikz(
      p = peak_bynetwork,
      file = "analyses/figures/peak-bynetwork.tex",
      width = 6,
      height = 4
    ),
    format = "file"
  ),
  tar_target(
    topo,
    make_topo(
      data_topo_gold = data_topo_gold,
      data_topo_gold_to_study = data_topo_gold_to_study,
      data_topo_study_to_study = data_topo_study_to_study,
      data_topo_sub_to_sub = data_topo_sub_to_sub
    ),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    fig_topo,
    make_tikz(p = topo, file = "analyses/figures/topo.tex", width = 3.25, height = 7),
    format = "file"
  ),
  tar_target(
    prop_effect_size,
    make_prop_effect_size(data_topo_gold = data_topo_gold),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    fig_prop_effect_size,
    make_tikz(
      p = prop_effect_size,
      file = "analyses/figures/prop-effect-size.tex",
      width = 4.5,
      height = 3
    ),
    format = "file"
  ),
  tar_target(
    topo_bynetwork,
    make_topo_bynetwork(pop_cor_region = pop_cor_region, data_topo_gold, at),
    packages = c("ggplot2")
  ),
  tar_target(
    fig_topo_bynetwork,
    make_tikz(
      p = topo_bynetwork,
      file = "analyses/figures/topo-bynetwork.tex",
      width = 6,
      height = 4
    ),
    format = "file"
  ),
  tar_target(
    model,
    make_model(
      data_model_gold_gold_to_study = data_model_gold_gold_to_study,
      data_model_study_to_study = data_model_study_to_study,
      data_model_sub_to_sub = data_model_sub_to_sub
    ),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    fig_model,
    make_tikz(p = model, file = "analyses/figures/model.tex", width = 3.25, height = 6),
    format = "file"
  ),
  tar_target(
    all_cog,
    make_all_cog(data_model_gold_gold_to_study = data_model_gold_gold_to_study),
    packages = c("ggplot2")
  ),
  tar_target(
    fig_all_cog,
    make_tikz(
      p = all_cog,
      file = "analyses/figures/all_cog.tex",
      width = 6,
      height = 8
    ),
    format = "file"
  ),
  tar_target(
    model_all_consistency,
    make_model_all_icc(
      data_model_study_to_study = data_model_study_to_study,
      type = "consistency"
    ),
    packages = c("ggplot2")
  ),
  tar_target(
    model_all_agreement,
    make_model_all_icc(
      data_model_study_to_study = data_model_study_to_study,
      type = "agreement"
    ),
    packages = c("ggplot2")
  ),
  tar_target(
    fig_model_all_agreement,
    make_tikz(
      p = model_all_agreement,
      file = "analyses/figures/model_all_agreement.tex",
      width = 6,
      height = 8
    ),
    format = "file"
  ),
  tar_target(
    fig_model_all_consistency,
    make_tikz(
      p = model_all_consistency,
      file = "analyses/figures/model_all_consistency.tex",
      width = 6,
      height = 8
    ),
    format = "file"
  )
)
