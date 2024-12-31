
library(dplyr)
library(readr)
library(tidyr)

avail <- targets::tar_read(test) |>
  unnest(avail) |>
  mutate(
    Subject = stringr::str_extract(avail, "[[:digit:]]{6}") |> as.integer()
  ) |>
  distinct(Subject)
  

twins <- readr::read_csv("data-raw/hcp/restricted.csv") |>
  semi_join(avail) |>
  select(Subject, ZygositySR, ZygosityGT, Family_ID, Mother_ID, Father_ID)


twins_s1200 <- readr::read_csv("data-raw/hcp/restricted.csv") |>
  select(Subject, ZygositySR, ZygosityGT, Family_ID, Mother_ID, Father_ID)


twins |>
  count(ZygositySR)

twins_s1200 |>
  count(ZygositySR)


twins_s1200 |>
  mutate(
    ZygositySR = case_match(
      ZygositySR,
      "NotMZ" ~ "DZ",
      .default = ZygositySR
    )
  ) |>
  count(ZygositySR, ZygosityGT)

