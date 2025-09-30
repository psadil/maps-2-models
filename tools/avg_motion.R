library(dplyr)

m <- readr::read_tsv("data-raw/motion.tsv") |> filter(!filtered)

test <- targets::tar_read(test) |>
  mutate(sub = as.numeric(sub), task = stringr::str_to_lower(Task))

test |>
  left_join(m, by = join_by(sub, task)) |>
  summarise(
    a = mean(loc, na.rm = TRUE),
    s = sd(loc, na.rm = TRUE)
  )

test_ukb <- targets::tar_read(test_ukb) |>
  mutate(sub = as.numeric(sub), task = stringr::str_to_lower(Task)) |>
  select(-UKB)


test_ukb |>
  left_join(select(m, -task), by = join_by(sub)) |>
  summarise(
    a = mean(loc, na.rm = TRUE),
    s = sd(loc, na.rm = TRUE)
  )
