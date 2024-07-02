import argparse
from pathlib import Path
import polars as pl
import typing
import logging

from sklearn.pipeline import make_pipeline
from sklearn import (
    feature_selection,
    preprocessing,
    metrics,
    linear_model,
)

from scipy import stats
from scipy.stats._resampling import PermutationTestResult

import numpy as np

from mlxtend.evaluate import permutation_test

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


def do_sample(
    clf,
    X_train: pl.DataFrame,
    y_train: pl.DataFrame,
    X_test: pl.DataFrame,
    y_test: pl.DataFrame,
) -> pl.DataFrame:
    model = clf.fit(X_train, y_train.to_numpy().squeeze())
    y_hat = model.predict(X_test)

    return y_test.with_columns(pl.Series("y_hat", y_hat))


def cor_rank(x, y) -> float:
    if (x[0] == x).all() or (y[0] == y).all():
        # If an input is constant, the correlation coefficient
        # is not defined.
        return np.NaN

    return stats.spearmanr(x, y).statistic


def test_cor(y, y_hat, seed: int | None = None) -> PermutationTestResult:
    def statistic(x):  # permute only `x`
        return stats.spearmanr(x, y).statistic

    return stats.permutation_test(
        (y_hat,), statistic, permutation_type="pairings", random_state=seed
    )


def cross_validate0(
    clf, X, y: pl.DataFrame, sub: pl.DataFrame, cv=None
) -> pl.DataFrame:
    model = clf.fit(X, y.to_series())
    y_hat = pl.Series("y_hat", model.predict(X))
    out = y.with_columns(y_hat).with_columns(
        sub=sub.to_series(),
        fold=0,
        statistic=cor_rank(y["g"], y_hat),
        r2=metrics.r2_score(y["g"], y_hat),
        mae=metrics.mean_absolute_error(y["g"], y_hat),
    )
    return out


def test_sample(
    d_test: pl.DataFrame,
    d_trainval: pl.DataFrame,
    clf,
    seed: int,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    out_ = cross_validate0(
        clf=clf,
        X=d_trainval.drop("sub", "g"),
        y=d_trainval.select("g"),
        sub=d_trainval.select("sub"),
    )

    # final train on full set
    model = clf.fit(d_trainval.drop("sub", "g"), d_trainval["g"])
    y_hat = model.predict(d_test.drop("sub", "g"))

    p = permutation_test(
        d_test["g"],
        y_hat,
        method="approximate",
        num_rounds=1000,
        func=lambda x, y: metrics.r2_score(x, y),
        seed=seed,
        paired=True,
    )

    result = test_cor(d_test["g"], y_hat, seed=seed)
    out = out_.with_columns(
        statistic_rep=result.statistic,
        pvalue_rep=result.pvalue,
        r2_rep_p=p,
        r2_rep=metrics.r2_score(d_test["g"], y_hat),
        mae_rep=metrics.median_absolute_error(d_test["g"], y_hat),
    )

    logging.info(f"{result.statistic=}")
    logging.info(f'r2: {metrics.r2_score(d_test["g"], y_hat)}')
    logging.info(f"{result.pvalue=}")
    out2 = d_test.select("sub", "g").with_columns(
        pl.Series(name="y_hat", values=y_hat), pl.col("g").cast(pl.Float32)
    )

    return out, out2


def main(
    x_in: Path,
    y_in: Path,
    out_dir: Path,
    n_studies: int,
    n_subs: list[int],
    tasks: list[TASK],
    measures: list[MEASURE],
) -> None:
    di = 64
    model = "RIDGE_CV"
    for measure in measures:
        for confound in [True, False]:
            logging.info(f"{measure=}")
            # https://github.com/KamalakerDadi/DiFuMo_analysis_scripts/blob/eeedd19b31e6e8859ba00fdbfb38b14e8fc88eea/3_3_Decoding_stimuli/ridge/run_pipeline_decoding_emotion.py#L60C24-L60C48
            clf = make_pipeline(
                feature_selection.VarianceThreshold(
                    0.01
                ),  # not all regions in fov
                preprocessing.RobustScaler(),
                linear_model.RidgeCV(alphas=np.logspace(-1.0, 4.0, 20)),
            )

            for task in tasks:
                logging.info(f"{task=}")
                y0 = pl.read_parquet(
                    y_in, columns=["Subject", measure]
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
                )

                d_test = d.sample(
                    n=int(d.shape[0] * 0.2),
                    shuffle=True,
                    seed=0,
                )

                d_trainval = d.join(d_test, on="sub", how="anti")

                out_, out2 = test_sample(
                    d_test=d_test,
                    d_trainval=d_trainval,
                    clf=clf,
                    seed=0,
                )

                parent = (
                    out_dir
                    / "out-perm-gold-cpm-sametest"
                    / f"measure={measure}"
                    / f"dimension={di}"
                    / f"confounds={confound}"
                    / f"model={model}"
                    / f"task={task}"
                )
                if not parent.exists():
                    parent.mkdir(parents=True)
                out_.write_parquet(parent / "part-0.parquet")
                parent = (
                    out_dir
                    / "out-perm-gold-cpm-preds-sametest"
                    / f"measure={measure}"
                    / f"dimension={di}"
                    / f"confounds={confound}"
                    / f"model={model}"
                    / f"task={task}"
                )
                if not parent.exists():
                    parent.mkdir(parents=True)
                out2.write_parquet(parent / "part-0.parquet")

                for n_sub in n_subs:
                    logging.info(f"{n_sub=}")
                    out = []
                    out2 = []
                    for study in range(n_studies):
                        logging.info(f"{study=}")

                        d_trainval = d.join(
                            d_test, on="sub", how="anti"
                        ).sample(
                            n=n_sub,
                            with_replacement=True,
                            shuffle=True,
                            seed=study,
                        )

                        out_, out2_ = test_sample(
                            d_test=d_test,
                            d_trainval=d_trainval,
                            clf=clf,
                            seed=study,
                        )

                        out.append(out_.with_columns(study=study))
                        out2.append(out2_.with_columns(study=study))

                    parent = (
                        out_dir
                        / "out-perm-cpm-sametest"
                        / f"measure={measure}"
                        / f"dimension={di}"
                        / f"confounds={confound}"
                        / f"model={model}"
                        / f"task={task}"
                        / f"n_sub={n_sub}"
                    )
                    if not parent.exists():
                        parent.mkdir(parents=True)
                    pl.concat(out).write_parquet(parent / "part-0.parquet")
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
                    pl.concat(out2).write_parquet(parent / "part-0.parquet")


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
    )

    # main(
    #     n_studies=args.n_studies,
    #     n_subs=args.n_subs,
    #     tasks=args.tasks,
    #     measures=args.measures,
    # )
