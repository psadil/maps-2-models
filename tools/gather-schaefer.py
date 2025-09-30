from pathlib import Path

import polars as pl

root = Path("../data-raw/cpm-schaefer2")
for pq in list(root.rglob("*parquet")):
    pl.read_parquet(pq).with_columns(
        pl.col("components_z").cast(pl.Float64)
    ).write_parquet(pq)


pl.scan_parquet(root).drop("components").with_columns(
    pl.col("type").cast(pl.Enum(categories=["MSMALL", "SURFACE", "VOL", "UKB"])),
    pl.col("task").cast(
        pl.Enum(
            categories=[
                "EMOTION",
                "GAMBLING",
                "LANGUAGE",
                "MOTOR",
                "RELATIONAL",
                "SOCIAL",
                "WM",
            ]
        )
    ),
).collect().pivot(
    values="components_z", on="i", index=["confounds", "sub", "task", "type"]
).write_parquet("../data-raw/cpm-schaefer.parquet")
