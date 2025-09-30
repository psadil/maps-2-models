library(dplyr)
library(tidyr)

targets::tar_load(c(roi_avg, roi_avg_ukb))

roi_avg <- roi_avg |> dplyr::filter(n_parcels == 400)
roi_avg_ukb <- roi_avg_ukb |> dplyr::filter(n_parcels == 400)

rois <- bind_rows(roi_avg, roi_avg_ukb)
rois |> arrow::write_parquet("data-raw/rois0.parquet")

# gold <- rois |>
#   summarise(
#     d = abs(mean(value) / sd(value)),
#     .by = c(type, Task, label)
#   ) |>
#   mutate(r = rank(d, ties.method = "first"), .by = c(type, Task)) |>
#   filter(r <= 10)

# rois |>
#   inner_join(gold, by = join_by(type, Task, label)) |>
#   select(sub, type, task = Task, value, r) |>
#   pivot_wider(names_from = r, names_sort = TRUE) |>
#   mutate(confounds = FALSE, sub = as.integer(sub)) |>
#   arrow::write_parquet("data-raw/rois-for-prediction.parquet")
