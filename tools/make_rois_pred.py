import polars as pl

roi_avg = pl.read_parquet("../data-raw/rois0.parquet")

gold = (
    roi_avg.group_by("type", "Task", "label")
    .agg(d=pl.col("value").mean().abs() / pl.col("value").std())
    .with_columns(r=pl.col("d").rank(method="ordinal").over("type", "Task"))
    .filter(pl.col("r") < 11)
)

roi_avg.join(gold, how="inner", on=["type", "Task", "label"]).rename(
    {"Task": "task"}
).select("sub", "task", "type", "value", "r").sort("r").pivot(
    on="r", index=["sub", "task", "type"]
).with_columns(
    pl.col("sub").cast(pl.Int64),
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
    confounds=False,
).write_parquet("../data-raw/rois-for-prediction.parquet")
