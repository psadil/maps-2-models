import polars as pl

u = pl.read_csv(
    "../data-raw/act_preds/data/unrestricted_psadil_2_20_2023_14_49_6.csv"
)
r = pl.read_csv(
    "../data-raw/act_preds/data/RESTRICTED_martin_5_16_2021_16_53_10.csv"
)

u.join(r, on="Subject").write_parquet("../data-raw/hcp.parquet")
