import argparse
import logging
import typing
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats
from scipy.stats._resampling import PermutationTestResult
from sklearn import (
    cross_decomposition,
    decomposition,
    feature_selection,
    linear_model,
    metrics,
    preprocessing,
)
from sklearn.pipeline import make_pipeline

logging.basicConfig(
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    level=logging.INFO,
)
ROOT = Path("/fastscratch/myscratch/pssadil")
# ROOT = Path("/Users/psadil/git/manuscripts/maps-to-models/meta/data-raw")
HCP_SUBS = ROOT / "hcp-subs.parquet"
UKB_SUBS = ROOT / "ukb-subs.parquet"
UKB_SMALL = ROOT / "ukb-sub-set.parquet"
HCP_Y = ROOT / "hcp.parquet"
UKB_Y = ROOT / "cognitive.parquet"

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
    "EMOTION", "GAMBLING", "LANGUAGE", "MOTOR", "RELATIONAL", "SOCIAL", "WM"
]

TYPE = typing.Literal["MSMALL", "SURFACE", "VOL", "UKB", "UKB_SMALL"]
MODEL = typing.Literal["RIDGE_CV", "PCR_RIDGE", "LASSO", "PCR_LASSO", "PLSR"]


def cor_rank(x, y) -> float:
    if (x[0] == x).all() or (y[0] == y).all():
        # If an input is constant, the correlation coefficient
        # is not defined.
        return np.nan

    return stats.spearmanr(x, y).statistic  # type: ignore


def test_cor(y, y_hat, seed: int | None = None) -> PermutationTestResult:
    def statistic(x):  # permute only `x`
        return stats.spearmanr(x, y).statistic  # type: ignore

    return stats.permutation_test(
        (y_hat,), statistic, permutation_type="pairings", rng=seed
    )


def test_r2(y, y_hat, seed: int | None = None) -> PermutationTestResult:
    def statistic(x):  # permute only `x`
        return metrics.r2_score(y, x)

    return stats.permutation_test(
        (y_hat,), statistic, permutation_type="pairings", rng=seed
    )


def test_sample(
    d_test: pl.DataFrame, d_trainval: pl.DataFrame, seed: int, m: MODEL = "RIDGE_CV"
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    # https://github.com/KamalakerDadi/DiFuMo_analysis_scripts/blob/eeedd19b31e6e8859ba00fdbfb38b14e8fc88eea/3_3_Decoding_stimuli/ridge/run_pipeline_decoding_emotion.py#L60C24-L60C48
    match m:
        case "RIDGE_CV":
            clf = make_pipeline(
                feature_selection.VarianceThreshold(0.01),  # not all regions in fov
                preprocessing.RobustScaler(),
                linear_model.RidgeCV(alphas=np.logspace(-1.0, 4.0, 20)),
            )
        case "LASSO":
            clf = make_pipeline(
                feature_selection.VarianceThreshold(0.01),  # not all regions in fov
                preprocessing.RobustScaler(),
                linear_model.LassoCV(),
            )
        case "PCR_RIDGE":
            clf = make_pipeline(
                feature_selection.VarianceThreshold(0.01),  # not all regions in fov
                preprocessing.RobustScaler(),
                decomposition.PCA(n_components=20),
                linear_model.RidgeCV(alphas=np.logspace(-1.0, 4.0, 20)),
            )
        case "PCR_LASSO":
            clf = make_pipeline(
                feature_selection.VarianceThreshold(0.01),  # not all regions in fov
                preprocessing.RobustScaler(),
                decomposition.PCA(n_components=20),
                linear_model.LassoCV(),
            )
        case "PLSR":
            clf = make_pipeline(
                feature_selection.VarianceThreshold(0.01),  # not all regions in fov
                cross_decomposition.PLSRegression(n_components=20),
            )
        case _:
            raise AssertionError("Uknown Model")

    model = clf.fit(d_trainval.drop("sub", "g"), d_trainval["g"])
    y_hat_trainvals = pl.Series("y_hat", model.predict(d_trainval.drop("sub", "g")))
    out_ = (
        d_trainval.select("g")
        .with_columns(y_hat_trainvals)
        .with_columns(
            sub=d_trainval["sub"],
            rank_cor=cor_rank(d_trainval["g"], y_hat_trainvals),
            r2=metrics.r2_score(d_trainval["g"], y_hat_trainvals),  # type: ignore
            mae=metrics.mean_absolute_error(d_trainval["g"], y_hat_trainvals),  # type: ignore
        )
    )

    y_hat = model.predict(d_test.drop("sub", "g"))

    result_r2 = test_r2(d_test["g"], y_hat, seed=seed)
    result = test_cor(d_test["g"], y_hat, seed=seed)
    out = out_.with_columns(
        statistic_rep=result.statistic,
        pvalue_rep=result.pvalue,
        r2_rep_p=result_r2.pvalue,
        r2_rep=result_r2.statistic,
        mae_rep=metrics.median_absolute_error(d_test["g"], y_hat),  # type: ignore
    )

    logging.info(f"cor: {result.statistic=}")
    logging.info(f"cor: {result.pvalue=}")
    logging.info(f"r2: {result_r2.statistic=}")
    logging.info(f"r2: {result_r2.pvalue=}")
    out2 = d_test.select("sub", "g").with_columns(
        pl.Series(name="y_hat", values=y_hat, dtype=pl.Float32),
        pl.col("g").cast(pl.Float32),
    )
    match m:
        case "RIDGE_CV" | "LASSO":
            coef: np.ndarray = model[2].coef_
        case "PCR_RIDGE" | "PCR_LASSO":
            coef: np.ndarray = model[3].coef_
        case "PLSR":
            coef: np.ndarray = model[1].coef_
        case _:
            raise AssertionError("Uknown Model")

    if not coef.dtype == np.float64:
        coef = np.asarray(coef, dtype=np.float64)
    if len(coef.shape) > 1:
        coef = coef.squeeze()

    return (out, out2, pl.DataFrame({"coef": coef}).with_row_index())


def main(
    x_in: Path,
    out_dir: Path,
    n_studies: int,
    n_subs: list[int],
    task: TASK,
    m: int,
    confound: bool,
    t: TYPE,
    model: MODEL = "RIDGE_CV",
    replacement: bool = True,
) -> None:
    logging.info(f"{task=}")
    logging.info(f"{t=}")
    logging.info(f"{confound=}")

    match t:
        case "MSMALL" | "SURFACE" | "VOL":
            filt = pl.scan_parquet(HCP_SUBS)
            measure = typing.get_args(MEASURE)[m]
            y0 = (
                pl.scan_parquet(HCP_Y)
                .select("sub", measure)
                .rename({measure: "g"})
                .drop_nulls()
            )
        case "UKB" | "UKB_SMALL":
            if t == "UKB":
                ukb_subs = UKB_SUBS
            elif t == "UKB_SMALL":
                ukb_subs = UKB_SMALL

            filt = pl.scan_parquet(ukb_subs).with_columns(
                task=pl.lit("EMOTION").cast(pl.Enum(categories=typing.get_args(TASK)))
            )
            ukb = pl.scan_parquet(UKB_Y).head(1).collect()
            measure = ukb.columns[m]
            y0 = (
                pl.scan_parquet(UKB_Y)
                .select("sub", measure)
                .rename({measure: "g"})
                .drop_nulls()
            )
        case _:
            raise ValueError("Unknown Type")

    logging.info(f"{measure=}")
    final_parent = (
        out_dir
        / "out-perm-cpm-preds-sametest-schaefer"
        / f"type={t}"
        / f"measure={measure}"
        / f"confounds={confound}"
        / f"replacement={replacement}"
        / f"model={model}"
        / f"task={task}"
        / f"n_sub={n_subs[-1]}"
    )
    if final_parent.exists():
        logging.info("skipping. outputs exist")
        return

    t2 = "UKB" if t == "UKB_SMALL" else t
    d = (
        pl.scan_parquet(x_in)
        .join(filt, how="semi", on=["sub", "task"])
        .filter(
            (pl.col("task") == task)
            & (pl.col("type") == t2)
            & (pl.col("confounds") == confound)
        )
        .drop("task", "type", "confounds")
        .join(y0, on="sub")
        .collect()
    )
    logging.info(f"N={d.height}")
    # drop columns that have all nulls (VOLS atlas may have different
    # numbers of regions in parcellation)
    d = d[[s.name for s in d if not (s.null_count() == d.height)]]
    d = d.drop_nans()

    d_test = d.sample(fraction=0.2, shuffle=True, seed=0)

    d_trainval = d.join(d_test, on="sub", how="anti")

    out_, out2_, features_ = test_sample(
        d_test=d_test, d_trainval=d_trainval, seed=0, m=model
    )

    parent = (
        out_dir
        / "out-perm-gold-cpm-sametest-schaefer"
        / f"type={t}"
        / f"measure={measure}"
        / f"confounds={confound}"
        / f"model={model}"
        / f"task={task}"
    )
    if not parent.exists():
        parent.mkdir(parents=True)
        out_.write_parquet(parent / "part-0.parquet")
    else:
        logging.info("gold already exists, skipping")
    parent = (
        out_dir
        / "out-perm-gold-cpm-preds-sametest-schaefer"
        / f"type={t}"
        / f"measure={measure}"
        / f"confounds={confound}"
        / f"model={model}"
        / f"task={task}"
    )
    if not parent.exists():
        parent.mkdir(parents=True)
        out2_.write_parquet(parent / "part-0.parquet")
    else:
        logging.info("gold2 already exists, skipping")

    parent = (
        out_dir
        / "out-perm-gold-cpm-preds-sametest-schaefer-features"
        / f"type={t}"
        / f"measure={measure}"
        / f"confounds={confound}"
        / f"model={model}"
        / f"task={task}"
    )
    if not parent.exists():
        parent.mkdir(parents=True)
        features_.write_parquet(parent / "part-0.parquet")
    else:
        logging.info("gold3 already exists, skipping")

    for n_sub in n_subs:
        logging.info(f"{n_sub=}")
        out = []
        out2 = []
        features = []
        for study in range(n_studies):
            logging.info(f"{study=}")

            out_, out2_, features_ = test_sample(
                d_test=d_test,
                d_trainval=d_trainval.sample(
                    n=n_sub, with_replacement=replacement, shuffle=True, seed=study
                ),
                seed=study,
                m=model,
            )

            out.append(out_.with_columns(study=study))
            out2.append(out2_.with_columns(study=study))
            features.append(features_.with_columns(study=study))

        parent = (
            out_dir
            / "out-perm-cpm-sametest-schaefer"
            / f"type={t}"
            / f"measure={measure}"
            / f"confounds={confound}"
            / f"replacement={replacement}"
            / f"model={model}"
            / f"task={task}"
            / f"n_sub={n_sub}"
        )
        if not parent.exists():
            parent.mkdir(parents=True)
        pl.concat(out).write_parquet(parent / "part-0.parquet")
        parent = (
            out_dir
            / "out-perm-cpm-preds-sametest-schaefer"
            / f"type={t}"
            / f"measure={measure}"
            / f"confounds={confound}"
            / f"replacement={replacement}"
            / f"model={model}"
            / f"task={task}"
            / f"n_sub={n_sub}"
        )
        if not parent.exists():
            parent.mkdir(parents=True)
        pl.concat(out2).write_parquet(parent / "part-0.parquet")
        parent = (
            out_dir
            / "out-perm-cpm-preds-sametest-schaefer-features"
            / f"type={t}"
            / f"measure={measure}"
            / f"confounds={confound}"
            / f"replacement={replacement}"
            / f"model={model}"
            / f"task={task}"
            / f"n_sub={n_sub}"
        )
        if not parent.exists():
            parent.mkdir(parents=True)
        pl.concat(features).write_parquet(parent / "part-0.parquet")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("out", type=Path)
    parser.add_argument("x", type=Path)
    parser.add_argument("measure", type=int)
    parser.add_argument("type", choices=typing.get_args(TYPE))
    parser.add_argument("task", choices=typing.get_args(TASK))
    parser.add_argument(
        "--confounds", action=argparse.BooleanOptionalAction, default=False
    )
    parser.add_argument("--model", default="RIDGE_CV", choices=typing.get_args(MODEL))

    parser.add_argument("--n-studies", type=int, default=100)
    parser.add_argument("--n-subs", nargs="+", type=int, default=[20, 40, 60, 80, 100])
    parser.add_argument(
        "--replacement", action=argparse.BooleanOptionalAction, default=True
    )
    args = parser.parse_args()

    main(
        x_in=args.x,
        out_dir=args.out,
        n_studies=args.n_studies,
        n_subs=args.n_subs,
        task=args.task,
        m=args.measure,
        confound=args.confounds,
        t=args.type,
        model=args.model,
        replacement=args.replacement,
    )
