import argparse
from pathlib import Path
import polars as pl
import typing
import logging

from sklearn.pipeline import make_pipeline
from sklearn import feature_selection, preprocessing, linear_model
from sklearn import model_selection


import numpy as np


logging.basicConfig(
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    level=logging.INFO,
)


MEASURE = typing.Literal[
    "Age_in_Yrs",
    "PicSeq_AgeAdj",
    "CardSort_AgeAdj",
    "Flanker_AgeAdj",
    "ReadEng_AgeAdj",
    "PicVocab_AgeAdj",
    "ProcSpeed_AgeAdj",
    "IWRD_TOT",
    "ListSort_AgeAdj",
    "CogTotalComp_AgeAdj",
    "CogCrystalComp_AgeAdj",
    "CogFluidComp_AgeAdj",
    "PMAT24_A_CR",
    "ASR_Extn_T",
    "ASR_Intn_T",
    "ASR_Attn_T",
    "NEOFAC_O",
    "NEOFAC_C",
    "NEOFAC_E",
    "NEOFAC_A",
    "NEOFAC_N",
    "DDisc_AUC_40K",
    "SCPT_SEN",
    "SCPT_SPEC",
    "VSPLOT_TC",
    "MMSE_Score",
    "PSQI_Score",
    "Endurance_Unadj",
    "GaitSpeed_Comp",
    "Dexterity_Unadj",
    "Strength_Unadj",
    "Odor_Unadj",
    "PainInterf_Tscore",
    "Taste_Unadj",
    "Mars_Final",
    "Emotion_Task_Face_Acc",
    "Language_Task_Math_Avg_Difficulty_Level",
    "Language_Task_Story_Avg_Difficulty_Level",
    "Social_Task_Perc_Random",
    "Social_Task_Perc_TOM",
    "WM_Task_Acc",
    "ER40_CR",
    "ER40FEAR",
    "ER40HAP",
    "ER40NOE",
    "ER40SAD",
    "AngAffect_Unadj",
    "AngHostil_Unadj",
    "AngAggr_Unadj",
    "FearAffect_Unadj",
    "FearSomat_Unadj",
    "Sadness_Unadj",
    "LifeSatisf_Unadj",
    "MeanPurp_Unadj",
    "PosAffect_Unadj",
    "Friendship_Unadj",
    "Loneliness_Unadj",
    "PercHostil_Unadj",
    "PercReject_Unadj",
    "EmotSupp_Unadj",
    "InstruSupp_Unadj",
    "PercStress_Unadj",
    "SelfEff_Unadj",
]

TASK = typing.Literal[
    "WM",
    "GAMBLING",
    "MOTOR",
    "LANGUAGE",
    "RELATIONAL",
    "SOCIAL",
    "EMOTION",
]


def test_sample(
    d_test: pl.DataFrame, d_trainval: pl.DataFrame
) -> pl.DataFrame:
    # https://github.com/KamalakerDadi/DiFuMo_analysis_scripts/blob/eeedd19b31e6e8859ba00fdbfb38b14e8fc88eea/3_3_Decoding_stimuli/ridge/run_pipeline_decoding_emotion.py#L60C24-L60C48
    # clf = make_pipeline(
    #     feature_selection.VarianceThreshold(0.01),  # not all regions in fov
    #     preprocessing.RobustScaler(),
    #     linear_model.RidgeCV(alphas=np.logspace(-1.0, 4.0, 20)),
    # )
    clf = make_pipeline(
        feature_selection.VarianceThreshold(0.01),  # not all regions in fov
        preprocessing.RobustScaler(),
        linear_model.LassoCV(cv=model_selection.KFold(n_splits=5)),
    )

    model = clf.fit(d_trainval.drop("g"), d_trainval["g"])
    y_hat = model.predict(d_test.drop("sub", "g"))

    out2 = d_test.select("sub", "g").with_columns(
        pl.Series(name="y_hat", values=y_hat), pl.col("g").cast(pl.Float32)
    )

    return out2


def main(
    x_in: Path,
    y_in: Path,
    out_dir: Path,
    n_studies: int,
    n_subs: list[int],
    tasks: list[TASK],
    measures: list[MEASURE],
    confounds: list[bool],
) -> None:
    di = 64
    model = "LASSO"
    outer_cv = model_selection.GroupShuffleSplit(
        n_splits=1, test_size=0.2, random_state=0
    )
    for measure in measures:
        for confound in confounds:
            logging.info(f"{measure=}")

            for task in tasks:
                logging.info(f"{task=}")
                y0 = pl.read_parquet(
                    y_in, columns=["Subject", measure, "Family_ID"]
                ).rename({measure: "g", "Subject": "sub"})

                d = (
                    pl.scan_parquet(x_in)
                    .filter(
                        (pl.col("task") == task)
                        & (pl.col("dimension") == di)
                        & (pl.col("confounds") == confound)
                    )
                    .drop("task", "dimension", "confounds")
                    .collect()
                    .join(y0, on="sub")
                    .drop_nulls()
                    .with_row_index()
                )

                for study, (_, test_index) in enumerate(
                    outer_cv.split(
                        d,  # type: ignore
                        groups=d.select("Family_ID").to_series(),
                    )
                ):
                    d_test = d.filter(pl.col("index").is_in(test_index)).drop(
                        "Family_ID", "index"
                    )
                    for n_sub in n_subs:
                        logging.info(f"{n_sub=}")
                        out2 = []
                        # for study, (train_index, test_index) in enumerate(
                        #     outer_cv.split(
                        #         d,  # type: ignore
                        #         groups=d.select("Family_ID").to_series(),
                        #     )
                        # ):
                        for study in range(n_studies):

                            logging.info(f"{study=}")
                            # d_test = d.filter(
                            #     pl.col("index").is_in(test_index)
                            # ).drop("Family_ID", "index")
                            d_trainval = (
                                d.join(d_test, on="sub", how="anti")
                                .sample(
                                    n=n_sub,
                                    with_replacement=True,
                                    shuffle=True,
                                    seed=study,
                                )
                                .drop("index", "Family_ID", "sub")
                            )

                            # d_trainval = (
                            #     d.filter(pl.col("index").is_in(train_index))
                            #     .drop("index", "Family_ID", "sub")
                            #     .sample(
                            #         n=n_sub,
                            #         with_replacement=True,
                            #         shuffle=True,
                            #         seed=study,
                            #     )
                            # )

                            out2_ = test_sample(
                                d_test=d_test, d_trainval=d_trainval
                            )

                            out2.append(out2_.with_columns(study=study))

                        parent = (
                            out_dir
                            / "out-perm-cpm-preds-sametest"
                            / f"measure={measure}"
                            / f"dimension={di}"
                            / f"confounds={confound}"
                            / f"model={model}"
                            / f"task={task}"
                            / f"n_sub={n_sub}"
                        )
                        if not parent.exists():
                            parent.mkdir(parents=True)
                        pl.concat(out2).write_parquet(
                            parent / "part-0.parquet"
                        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("out", type=Path)
    parser.add_argument("x", type=Path)
    parser.add_argument("y", type=Path)
    parser.add_argument("--n-studies", type=int, default=100)
    parser.add_argument(
        "--n-subs", nargs="+", type=int, default=[20, 40, 60, 80, 100]
    )
    parser.add_argument(
        "--tasks",
        nargs="+",
        choices=typing.get_args(TASK),
        default=typing.get_args(TASK),
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--measures",
        nargs="+",
        choices=typing.get_args(MEASURE),
        default=typing.get_args(MEASURE),
    )
    group.add_argument("--measure", type=int)
    group.add_argument(
        "--confounds",
        nargs="+",
        choices=["True", "False"],
        default=["True", "False"],
    )
    args = parser.parse_args()

    if args.measure:
        measures = [typing.get_args(MEASURE)[args.measure]]
    else:
        measures = args.measures

    main(
        x_in=args.x,
        y_in=args.y,
        out_dir=args.out,
        n_studies=args.n_studies,
        n_subs=args.n_subs,
        tasks=args.tasks,
        measures=measures,
        confounds=[c == "True" for c in args.confounds],
    )
