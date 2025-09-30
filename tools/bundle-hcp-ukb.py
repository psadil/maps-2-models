import polars as pl

u = pl.read_csv("../data-raw/hcp/unrestricted.csv")
r = pl.read_csv("../data-raw/hcp/restricted.csv")

u.join(r, on="Subject").rename({"Subject": "sub"}).write_parquet(
    "../data-raw/hcp.parquet"
)

ukb = pl.read_csv(
    "../data-raw/cognitive.tsv", separator="\t", infer_schema_length=None
).rename({"f.eid": "sub"})
ukb = ukb[[s.name for s in ukb if not (s.null_count() == ukb.height)]]

ukb.select("sub", pl.selectors.matches(r"f\.\d+\.2\.\d+")).write_parquet(
    "../data-raw/cognitive.parquet"
)
