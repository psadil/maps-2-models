import polars as pl

u = pl.read_csv("../data-raw/hcp/unrestricted.csv")
r = pl.read_csv("../data-raw/hcp/restricted.csv")

u.join(r, on="Subject").write_parquet("../data-raw/hcp.parquet")
