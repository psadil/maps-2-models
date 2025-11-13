library(targets)
library(tarchetypes)
library(rlang)

source(here::here("R", "spatial.R"))
source(here::here("R", "tfce.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "hcp.R"))
source(here::here("R", "ptfce.R"))
source(here::here("R", "figures.R"))
source(here::here("R", "model.R"))
source(here::here("R", "manuscript.R"))
source(here::here("R", "roi.R"))
source(here::here("R", "topo.R"))
source(here::here("R", "ukb.R"))
source(here::here("R", "cifti.R"))
source(here::here("R", "peaks.R"))

Sys.setenv(
  NIIDIR = here::here("data-raw", "hcp-niis-ptfce"),
  PALMDIR = here::here("data-raw", "palm"),
  HCPPARQUET = "/Users/psadil/Library/CloudStorage/OneDrive-JohnsHopkins/data/hcp-to-parquet/data/out",
  VOLGLM = "/dcl01/smart/data/psadil/meta/data-raw/glm_manual",
  UKBMNI = "/dcl01/smart/data/psadil/meta/data-raw/ukb_mni"
) # explicitly avoiding tracking this

# controller_small <- crew::crew_controller_local(
#   name = "small",
#   workers = 1
# )

controller <- crew::crew_controller_local(
  name = "large",
  workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
)

controller_small <- crew::crew_controller_local(
  name = "small",
  workers = 1
)

# controller <- crew.cluster::crew_controller_slurm(
#   name = "slurm",
#   workers = 5,
#   tasks_max = 10,
#   options_cluster = crew.cluster::crew_options_slurm(
#     script_lines = "source ~/.bashrc; mamba activate meta; ml matlab; ml fsl; export PATH=/dcl01/smart/data/psadil/meta/tools:/users/pssadil/workbench/bin_rh_linux64:${PATH}:/users/pssadil/git/PALM",
#     memory_gigabytes_required = 2,
#     cpus_per_task = 1,
#     time_minutes = 10,
#     command_submit = '/usr/bin/sbatch --constraint="intel"'
#   )
# )

targets::tar_option_set(
  format = "qs",
  storage = "worker",
  packages = c("oro.nifti"),
  controller = crew::crew_controller_group(controller, controller_small),
  workspace_on_error = FALSE
)


list(
  tar_group_by(
    contrasts,
    get_hcp_contrasts(),
    Task,
    CopeNumber,
    format = "parquet",
    deployment = "main"
  ),
  tar_target(
    test,
    get_hcp_copes(contrasts),
    format = "parquet",
    deployment = "main"
  ),
  tar_target(
    test_ukb,
    get_ukb_copes("data-raw/ukb_copes"),
    format = "parquet",
    deployment = "main"
  ),
  tar_target(
    test_all,
    get_hcp_copes(contrasts, matching_only = FALSE),
    format = "parquet",
    deployment = "main"
  ),
  tar_target(n_sub, c(20, 40, 60, 80, 100), deployment = "main"),
  tar_target(iter, seq_len(100), deployment = "main"),
  tar_group_by(
    hcp_samples,
    sample_hcp(test = test, n_iter = 100, n_subs = c(20, 40, 60, 80, 100)),
    Task,
    CopeNumber,
    type,
    iter,
    n_sub,
    format = "parquet"
  ),
  tar_group_by(
    ukb_samples,
    sample_ukb(
      test = test_ukb,
      n_iter = 100,
      n_subs = c(20, 40, 60, 80, 100, 1000, 10000)
    ),
    Task,
    CopeNumber,
    type,
    iter,
    n_sub,
    format = "parquet"
  ),
  tar_group_by(
    hcp_samples_all,
    sample_hcp(
      test = test_all,
      n_iter = 100,
      n_subs = c(20, 40, 60, 80, 100),
      types = c("MSMALL")
    ),
    Task,
    CopeNumber,
    type,
    iter,
    n_sub,
    format = "parquet"
  ),
  tar_target(
    roi_avg,
    avg_roi(test, n_parcels),
    pattern = cross(test, n_parcels),
    format = "parquet"
  ),
  tar_target(rois, test_roi_wrapper(roi_avg, hcp_samples), format = "parquet"),
  tar_target(rois_pop, test_roi_pop(roi_avg), format = "parquet"),
  tar_target(
    roi_avg_ukb,
    avg_roi_ukb(test_ukb, n_parcels),
    pattern = cross(test_ukb, n_parcels),
    format = "parquet"
  ),
  tar_target(
    rois_ukb,
    test_roi_wrapper(roi_avg_ukb, ukb_samples),
    format = "parquet"
  ),
  tar_target(rois_pop_ukb, test_roi_pop(roi_avg_ukb), format = "parquet"),
  # strategy for palm: create commands, then run with slurm array
  tar_target(
    tfce_ukb,
    get_tfce_cmd(ukb_samples, test_ukb),
    format = "parquet",
    deployment = "main"
  ),
  tar_target(
    tfce,
    dplyr::bind_rows(get_tfce_cmd(hcp_samples, test), tfce_ukb),
    format = "parquet",
    deployment = "main"
  ),
  # this is back in R/targets (assumes above have been created)
  # (palm does not need to be run for tfce_pop{_ukb}, but we do need target)
  tar_target(tfce_pop_ukb, get_tfce_pop(test_ukb), format = "parquet"),
  tar_target(
    tfce_pop,
    dplyr::bind_rows(get_tfce_pop(test), tfce_pop_ukb),
    format = "parquet"
  ),
  tar_target(threshold, c(0, -log(0.05))),
  tar_target(
    study_peaks,
    get_study_peaks_cifti_rows(tfce, max_n_sub = 40, threshold = threshold),
    format = "parquet",
    pattern = map(threshold)
  ),
  tar_target(
    gold_peaks_cifti,
    get_gold_peaks_cifti(tfce_pop),
    format = "parquet",
    pattern = map(tfce_pop)
  ),
  tarchetypes::tar_group_by(
    study_peaks_grouped,
    study_peaks,
    type,
    Task,
    CopeNumber,
    n_sub,
    iter,
    threshold,
    format = "parquet",
    deployment = "main"
  ),
  tar_target(
    study_to_gold_distances,
    get_cifti_augmented2_rows(study_peaks_grouped, gold_peaks_cifti),
    format = "parquet",
    pattern = map(study_peaks_grouped)
  ),
  tar_target(
    glm,
    get_glm(test, hcp_samples),
    format = "parquet",
    pattern = map(hcp_samples)
  ),
  tar_target(
    glm_pop,
    get_glm_pop(test, contrasts),
    format = "parquet",
    pattern = map(contrasts)
  ),
  tar_target(
    glm_ukb,
    get_glm(
      dplyr::mutate(
        test_ukb,
        UKB = stringr::str_replace(
          UKB,
          "/fastscratch/myscratch/pssadil/ukb_mni/derivatives/",
          "data-raw/ukb_mni/"
        ) |>
          stringr::str_remove("_task-EMOTION_space-MNI152_contrast-5_cope") |>
          stringr::str_remove("sub-")
      ),
      ukb_samples
    ),
    format = "parquet",
    pattern = map(ukb_samples)
  ),
  tar_target(
    ukb_gray,
    prep_ukb_pop(test_ukb),
    format = "parquet",
    pattern = map(test_ukb)
  ),
  # ukb_glm_pop.parquet make with tools/run_ukb_pop_glm.R
  tar_target(
    glm_pop_ukb_file,
    "data-raw/ukb_glm_pop.parquet",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    glm_pop_ukb,
    get_glm_pop_ukb(glm_pop_ukb_file),
    format = "parquet"
  ),
  tar_target(glm2, dplyr::bind_rows(glm, glm_ukb), format = "parquet"),
  tar_target(
    glm_pop2,
    dplyr::bind_rows(glm_pop, glm_pop_ukb),
    format = "parquet"
  ),
  tar_target(
    data_topo_gold_to_study,
    make_data_topo_gold_to_study2(glm2, glm_pop2),
    format = "parquet"
  ),
  tar_target(
    data_topo_gold_to_study_bynetwork,
    make_data_topo_gold_to_study_bynetwork(
      glm2 = glm2,
      glm_pop2 = glm_pop2,
      at = at
    ),
    format = "parquet"
  ),
  tarchetypes::tar_group_by(
    data_topo_study_to_study0,
    glm2,
    Task,
    type,
    n_sub,
    format = "parquet"
  ),
  tar_target(
    data_topo_study_to_study,
    make_data_topo_study_to_study(data_topo_study_to_study0),
    pattern = map(data_topo_study_to_study0),
    format = "parquet"
  ),
  tar_target(at, make_atlas_full()),
  tar_target(n_parcels, c(200, 400, 800)),
  tar_target(
    ptfce,
    do_ptfce2(hcp_samples = hcp_samples, test = test, enhance = TRUE),
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
    get_ptfce_maxes(row = ptfce, do_fwe_correction = do_fwe_correction),
    pattern = cross(map(ptfce), do_fwe_correction),
    format = "parquet"
  ),
  tar_target(
    roi_avg_all,
    avg_roi(test_all, n_parcels),
    pattern = cross(test_all, n_parcels),
    format = "parquet"
  ),
  tar_target(
    rois_all,
    test_roi(dplyr::filter(roi_avg_all, n_parcels == 400), hcp_samples_all),
    pattern = map(hcp_samples_all),
    format = "parquet"
  ),
  tar_target(rois_pop_all, test_roi_pop(roi_avg_all), format = "parquet"),
  tar_target(
    data_peak_study_to_study,
    make_data_peak_study_to_study(study_to_gold_distances, at = at),
    format = "parquet"
  ),
  tar_target(
    data_roi_study_to_gold,
    make_data_roi_study_to_gold(
      gold_tested = dplyr::bind_rows(rois_pop, rois_pop_ukb),
      rois_tested = dplyr::bind_rows(rois, rois_ukb)
    ),
    format = "parquet"
  ),
  tar_target(
    data_roi_study_to_gold2,
    make_data_roi_study_to_gold2(
      gold_tested = dplyr::bind_rows(rois_pop, rois_pop_ukb),
      rois_tested = dplyr::bind_rows(rois, rois_ukb)
    ),
    format = "parquet"
  ),
  tar_target(
    data_roi_study_to_study,
    make_data_roi_study_to_study(
      rois_tested = dplyr::bind_rows(rois, rois_ukb)
    ),
    format = "parquet"
  ),
  tar_target(
    data_roi_study_to_study2,
    make_data_roi_study_to_study2(
      rois_tested = dplyr::bind_rows(rois, rois_ukb),
      n_workers = 8
    ),
    format = "parquet",
    resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "small")
    )
  ),
  tar_target(
    data_peak_study_to_gold,
    make_data_peak_study_to_gold(
      study_to_gold_distances,
      at = at,
      gold_tested = dplyr::bind_rows(rois_pop, rois_pop_ukb)
    ),
    format = "parquet"
  ),
  tar_target(
    data_model_gold_gold_to_study,
    make_data_model_gold_gold_to_study(
      dataset_gold = here::here("data-raw/out-perm-gold-cpm-sametest-schaefer"),
      dataset = here::here("data-raw/out-perm-cpm-sametest-schaefer"),
      measures = measures
    ),
    format = "parquet"
  ),
  tar_target(
    data_model_gold_gold_to_study2,
    make_data_model_gold_gold_to_study2(
      dataset_gold = here::here("data-raw/out-perm-gold-cpm-sametest-schaefer"),
      dataset = here::here("data-raw/out-perm-cpm-sametest-schaefer"),
      measures = measures
    ),
    format = "parquet"
  ),
  tar_target(
    data_model_gold_gold_to_study3,
    make_data_model_gold_gold_to_study3(
      dataset_gold = here::here(
        "data-raw/out-perm-gold-cpm-preds-sametest-schaefer-features"
      ),
      dataset = here::here(
        "data-raw/out-perm-cpm-preds-sametest-schaefer-features"
      ),
      measures = measures
    ),
    format = "parquet"
  ),
  tar_target(
    data_model_study_to_study,
    make_data_model_study_to_study(
      dataset = here::here("data-raw", "out-perm-cpm-preds-sametest-schaefer"),
      measures = measures
    ),
    format = "parquet"
  ),
  tar_target(
    data_model_study_to_study3,
    make_data_model_study_to_study3(
      dataset = here::here(
        "data-raw",
        "out-perm-cpm-preds-sametest-schaefer-features"
      ),
      measures = measures
    ),
    format = "parquet"
  ),
  tar_target(
    fig_roi,
    make_tikz(
      p = make_roi(data_roi_study_to_gold, data_roi_study_to_study),
      file = "analyses/figures/roi.tex",
      width = 6,
      height = 5
    ),
    format = "file",
    packages = c("patchwork")
  ),
  tar_target(
    fig_roi2,
    make_tikz(
      p = make_roi2(data_roi_study_to_gold2, data_roi_study_to_study2),
      file = "analyses/figures/roi2.tex",
      width = 7,
      height = 7.5
    ),
    format = "file",
    packages = c("patchwork")
  ),
  tar_target(
    fig_prop_active_most_active_roi,
    make_tikz(
      p = make_prop_active_most_active_roi(
        data_roi_study_to_gold2 = data_roi_study_to_gold2
      ),
      file = "analyses/figures/prop-active-most-active-roi.tex",
      width = 6,
      height = 8
    ),
    format = "file",
    packages = c("patchwork")
  ),
  tar_target(
    peaks_reliability,
    make_peaks_reliability(
      data_peak_study_to_study = data_peak_study_to_study,
      threshold = "reg",
      nrow_subfig = 2
    ),
    packages = c("patchwork")
  ),
  tar_target(
    peaks_validity,
    make_peaks_validity(
      data_peak_study_to_gold = data_peak_study_to_gold,
      threshold = "reg"
    ),
    packages = c("patchwork")
  ),
  tar_target(
    fig_peaks_validity,
    make_tikz(
      p = peaks_validity,
      file = "analyses/figures/peaks_validity.tex",
      width = 6,
      height = 8
    ),
    format = "file"
  ),
  tar_target(
    fig_peaks_reliability,
    make_tikz(
      p = peaks_reliability,
      file = "analyses/figures/peaks_reliability.tex",
      width = 7,
      height = 9
    ),
    format = "file"
  ),
  tar_target(
    topo,
    make_topo(
      data_topo_gold_to_study = data_topo_gold_to_study,
      data_topo_study_to_study = data_topo_study_to_study
    ),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    fig_topo,
    make_tikz(
      p = topo,
      file = "analyses/figures/topo.tex",
      width = 4.5,
      height = 6
    ),
    format = "file"
  ),
  tar_target(
    fig_prop_effect_size,
    make_tikz(
      p = make_prop_effect_size(glm_pop2),
      file = "analyses/figures/prop-effect-size.tex",
      width = 5,
      height = 5
    ),
    format = "file"
  ),
  tar_target(
    model,
    make_model(
      data_model_gold_gold_to_study = data_model_gold_gold_to_study,
      data_model_study_to_study = data_model_study_to_study
    ),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    fig_model,
    make_tikz(
      p = model,
      file = "analyses/figures/model.tex",
      width = 4.5,
      height = 6
    ),
    format = "file"
  ),
  tar_target(
    fig_model2,
    make_tikz(
      p = make_model2(
        data_model_gold_gold_to_study2 = data_model_gold_gold_to_study2
      ),
      file = "analyses/figures/model2.tex",
      width = 8,
      height = 6
    ),
    format = "file",
    packages = "patchwork"
  ),
  tar_target(
    fig_model2_neg,
    make_tikz(
      p = make_model2_neg(
        data_model_gold_gold_to_study2 = data_model_gold_gold_to_study2
      ),
      file = "analyses/figures/model2_neg.tex",
      width = 8,
      height = 6
    ),
    format = "file",
    packages = "patchwork"
  ),
  tar_target(
    fig_model3,
    make_tikz(
      p = make_model3(
        data_model_gold_gold_to_study3 = data_model_gold_gold_to_study3,
        data_model_study_to_study3 = data_model_study_to_study3
      ),
      file = "analyses/figures/model3.tex",
      width = 6,
      height = 8
    ),
    format = "file",
    packages = "patchwork"
  ),
  tar_target(
    fig_model_all_cog,
    make_tikz(
      p = make_all_cog(
        data_model_gold_gold_to_study = data_model_gold_gold_to_study,
        data_model_study_to_study = data_model_study_to_study
      ),
      file = "analyses/figures/model_all_cog.tex",
      width = 6,
      height = 8
    ),
    format = "file",
    packages = "patchwork"
  ),
  tar_target(
    model_model,
    make_model_model(
      data_model_gold_gold_to_study,
      data_model_gold_gold_to_study2,
      data_model_gold_gold_to_study3,
      data_model_study_to_study,
      data_model_study_to_study3
    ),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    model_model_ukb,
    make_model_model_ukb(
      data_model_gold_gold_to_study,
      data_model_gold_gold_to_study2,
      data_model_gold_gold_to_study3,
      data_model_study_to_study,
      data_model_study_to_study3
    ),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    model_sigmas,
    make_model_sigmas(data_model_study_to_study),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    model_sigmas_ukb,
    make_model_sigmas_ukb(data_model_study_to_study),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    model_sigmas3_ukb,
    make_model_sigmas3_ukb(data_model_study_to_study),
    packages = c("ggplot2", "patchwork")
  ),
  tar_target(
    peak_table_count,
    make_table_of_studies_without_peaks(
      study_peaks_grouped,
      "analyses/figures/peak_table_count.tex"
    ),
    format = "file"
  ),
  tar_target(
    fig_peaks_reliability_unthresholded,
    make_tikz(
      p = make_peaks_reliability(
        data_peak_study_to_study = data_peak_study_to_study,
        threshold = "unthresholded",
        nrow_subfig = 2,
        base_size = 8
      ),
      file = glue::glue(
        "analyses/figures/peaks_reliability_unthresholded.tex"
      ),
      width = 6,
      height = 8
    ),
    format = "file",
    packages = c("patchwork")
  ),
  tar_target(
    peaks_validity_unthresholded,
    make_peaks_validity(
      data_peak_study_to_gold = data_peak_study_to_gold,
      threshold = "unthresholded"
    ),
    packages = c("patchwork")
  ),
  tar_target(
    fig_peaks_validity_unthresholded,
    make_tikz(
      p = peaks_validity_unthresholded,
      file = "analyses/figures/peaks_validity_unthresholded.tex",
      width = 6,
      height = 8
    ),
    format = "file"
  ),
  tar_target(
    data_modelroi_gold_gold_to_study,
    make_data_model_gold_gold_to_study(
      dataset_gold = here::here(
        "data-raw/out-rois/out-perm-gold-cpm-sametest-schaefer"
      ),
      dataset = here::here("data-raw/out-rois/out-perm-cpm-sametest-schaefer"),
      measures = measures
    ),
    format = "parquet"
  ),
  tar_target(
    fig_modelroi,
    make_tikz(
      p = make_modelroi(data_modelroi_gold_gold_to_study),
      file = "analyses/figures/modelroi.tex",
      width = 4,
      height = 6
    ),
    format = "file"
  ),
  tar_target(
    data_model_gold_gold_to_study_r2,
    make_data_model_gold_gold_to_study_r2(
      dataset_gold = here::here(
        "data-raw/out-perm-gold-cpm-sametest-schaefer"
      ),
      dataset = here::here("data-raw/out-perm-cpm-sametest-schaefer"),
      measures = measures
    ),
    format = "parquet"
  ),
  tar_target(
    fig_model_r2,
    make_tikz(
      p = make_model_r2(data_model_gold_gold_to_study_r2),
      file = "analyses/figures/model_r2.tex",
      width = 6,
      height = 3
    ),
    format = "file"
  ),
  tar_target(measures, get_measures()),
  tar_target(
    fig_peaks_by_fwe,
    make_tikz(
      p = make_peaks_by_fwe(gold_peaks, maxes),
      file = "analyses/figures/peaks_by_fwe.tex",
      width = 8,
      height = 6
    ),
    format = "file",
    packages = c("patchwork")
  ),
  tar_target(
    fig_ecdf_peak_reliability,
    make_tikz(
      p = make_ecdf_peak_reliability(data_peak_study_to_study),
      file = "analyses/figures/ecdf_peak_reliability.tex",
      width = 6,
      height = 6
    ),
    format = "file",
    packages = c("patchwork")
  ),
  tar_target(
    top_ten_regions,
    write_regions(
      rois_pop,
      rois_pop_ukb,
      dst = "analyses/tables/top_ten_regions.tsv"
    ),
    format = "file"
  ),
  tar_target(
    performance_file,
    write_model_performance(
      data_model_gold_gold_to_study,
      data_model_gold_gold_to_study_r2,
      dst = "analyses/tables/mm_scores.tsv"
    ),
    format = "file"
  ),
  tar_target(
    top_ten_peaks,
    write_peaks(
      gold_tested = dplyr::bind_rows(rois_pop, rois_pop_ukb),
      dst = "analyses/tables/top_ten_peaks.tsv"
    ),
    format = "file"
  ),
  tar_target(
    fig_peak_bysize,
    make_tikz(
      p = make_peak_bysize(study_to_gold_distances, glm_pop2),
      file = "analyses/figures/peak-bysize.tex",
      width = 5,
      height = 5
    ),
    packages = "patchwork"
  ),
  tar_target(
    fig_peak_bynetwork,
    make_tikz(
      p = make_peak_bynetwork(study_to_gold_distances, at, glm_pop2),
      file = "analyses/figures/peak-bynetwork.tex",
      width = 5,
      height = 5
    ),
    packages = "patchwork"
  ),
  tar_target(
    fig_topo_bynetwork,
    make_tikz(
      p = make_topo_bynetwork(
        data_topo_gold_to_study_bynetwork,
        glm_pop2,
        at = at
      ),
      file = "analyses/figures/topo-bynetwork.tex",
      width = 5,
      height = 6
    ),
    packages = "patchwork"
  ),
  tar_target(
    peak_avg_bysize,
    make_peak_avg_bysize(study_to_gold_distances, glm_pop2)
  )
)
