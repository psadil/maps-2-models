import polars as pl
import pyarrow.dataset as ds

dset = ds.dataset(
    "/Users/psadil/git/manuscripts/maps-to-models/meta/data-raw/cpm-difumo2",
    format="parquet",
    partition_base_dir="/Users/psadil/git/manuscripts/maps-to-models/meta/data-raw/cpm-difumo2",
    partitioning="hive",
)

pl.scan_pyarrow_dataset(dset).collect().with_columns(
    pl.col("sub").cast(pl.Int64),
    confounds=pl.col("confounds").str.contains("True"),
).drop("components").pivot(
    values="components_z",
    on="i",
    index=["dimension", "confounds", "sub", "task"],
).write_parquet(
    "../data-raw/cpm-difumo.parquet"
)
