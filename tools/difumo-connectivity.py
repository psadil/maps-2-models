import argparse
from pathlib import Path
import re

from nilearn import datasets, connectome, image
from nilearn.maskers import NiftiMapsMasker
from sklearn import covariance
import polars as pl
import nibabel as nb
import pandas as pd
import numpy as np

OUTDIR = Path("/fastscratch") / "myscratch" / "pssadil" / "cpm-difumo2"

# create masker to extract functional data within atlas parcels
connectome_measure = connectome.ConnectivityMeasure(
    kind="correlation",
    cov_estimator=covariance.LedoitWolf(),
    vectorize=True,
    discard_diagonal=True,
)


def make_confounds(file: Path) -> np.ndarray:
    confounds = pl.from_dataframe(
        pd.read_csv(
            file.with_name("Movement_Regressors.txt"),
            sep=r"\s+",
            header=0,
            names=[
                "trans_x",
                "trans_y",
                "trans_z",
                "rot_x",
                "rot_y",
                "rot_z",
                "trans_x_derivative1",
                "trans_y_derivative1",
                "trans_z_derivative1",
                "rot_x_derivative1",
                "rot_y_derivative1",
                "rot_z_derivative1",
            ],
        )
    ).with_columns(
        trans_x_power2=pl.col("trans_x").pow(2),
        trans_y_power2=pl.col("trans_y").pow(2),
        trans_z_power2=pl.col("trans_z").pow(2),
        rot_x_power2=pl.col("rot_x").pow(2),
        rot_y_power2=pl.col("rot_y").pow(2),
        rot_z_power2=pl.col("rot_z").pow(2),
        trans_x_derivative1_power2=pl.col("trans_x_derivative1").pow(2),
        trans_y_derivative1_power2=pl.col("trans_y_derivative1").pow(2),
        trans_z_derivative1_power2=pl.col("trans_z_derivative1").pow(2),
        rot_x_derivative1_power2=pl.col("rot_x_derivative1").pow(2),
        rot_y_derivative1_power2=pl.col("rot_y_derivative1").pow(2),
        rot_z_derivative1_power2=pl.col("rot_z_derivative1").pow(2),
    )
    return np.vstack(
        [np.zeros(shape=(1, confounds.shape[1])), confounds.to_numpy()]
    )


def get_connectivity(img: Path, mask: Path, confounds: bool = True):
    for dim in [64]:
        difumo = datasets.fetch_atlas_difumo(
            dimension=dim, resolution_mm=2, legacy_format=False
        )
        masker = NiftiMapsMasker(maps_img=difumo.maps, standardize=False)
        print(f"{img=}")
        sub = re.findall(r"(?<=disk\d/)\d+", str(img))[0]
        task = re.findall(r"(?<=tfMRI_)[A-Z]+", str(img))[0]

        outfile = (
            OUTDIR
            / f"dimension={dim}"
            / f"confounds={confounds}"
            / f"sub={sub}"
            / f"task={task}"
            / "part-0.parquet"
        )
        if not (parent := outfile.parent).exists():
            parent.mkdir(parents=True)

        nii = nb.loadsave.load(img)
        if confounds:
            cnf_tbl = make_confounds(img)
        else:
            cnf_tbl = None
        cleaned = image.clean_img(
            imgs=img,
            high_pass=0.01,
            low_pass=0.1,
            detrend=True,
            standardize=True,
            confounds=cnf_tbl,
            t_r=nii.header.get("pixdim")[4],
            mask_img=mask,
        )
        if (rl_img := Path(str(img).replace("LR", "RL"))).exists():
            if confounds:
                cnf_tbl = make_confounds(rl_img)
            else:
                cnf_tbl = None
            rl_cleaned = image.clean_img(
                imgs=rl_img,
                high_pass=0.01,
                low_pass=0.1,
                detrend=True,
                standardize=True,
                confounds=cnf_tbl,
                t_r=nii.header.get("pixdim")[4],
                mask_img=mask,
            )
            cleaned = image.concat_imgs([cleaned, rl_cleaned])

        del nii
        x = masker.fit_transform(cleaned)
        del cleaned

        correlation_matrices = connectome_measure.fit_transform([x])
        pl.DataFrame({"components": correlation_matrices[0, :]}).with_columns(
            components_z=pl.col("components").arctanh(),
            i=pl.arange(0, correlation_matrices.shape[1]),
        ).write_parquet(outfile, statistics=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("i", type=int)  # total: 1225
    parser.add_argument(
        "--in-dir",
        type=Path,
        default=Path("/dcs04/legacy-dcs01-oasis/hpc/"),
    )
    parser.add_argument(
        "--mask",
        type=Path,
        default=Path(
            "/fastscratch/myscratch/pssadil/MNI152_T1_2mm_brain_mask_dil.nii.gz"
        ),
    )
    parser.add_argument(
        "--confounds", default=True, action=argparse.BooleanOptionalAction
    )

    args = parser.parse_args()
    to_process: list[Path] = sorted(
        list(
            args.in_dir.glob(
                "disk*/*/MNINonLinear/Results/*LR/tfMRI*LR.nii.gz"
            )
        )
    )

    get_connectivity(
        to_process[args.i], mask=args.mask, confounds=args.confounds
    )
