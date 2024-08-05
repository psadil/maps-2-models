import argparse
from pathlib import Path
import re

import numpy as np
import nibabel as nb
from nilearn import image

# from nibabel import affines
from nibabel import processing

import polars as pl

import prefect

# from prefect.task_runners import SequentialTaskRunner
from prefect_dask import DaskTaskRunner


FILES = {
    "cope1.nii.gz": np.float64,
    "mask.nii.gz": np.bool_,
    "pe1.nii.gz": np.float64,
    "tdof_t1.nii.gz": np.uint16,
    "varcope1.nii.gz": np.float64,
    "zflame1lowertstat1.nii.gz": np.float64,
    "zstat1.nii.gz": np.float64,
    "mean_random_effects_var1.nii.gz": np.float64,
    "tstat1.nii.gz": np.float64,
    "zflame1uppertstat1.nii.gz": np.float64,
}

ANAT = {"aparc.a2009s+aseg.nii.gz", "aparc+aseg.nii.gz"}


def _remove_niigz(f: Path) -> str:
    return f.name.removesuffix(".gz").removesuffix(".nii")


def to_polars0(f_stem: str, n: nb.Nifti1Image, dtype: np.dtype) -> pl.DataFrame:
    values = n.get_fdata().ravel()
    (i, j, k) = np.unravel_index(np.arange(len(values)), shape=n.shape)
    # xyz = affines.apply_affine(n.affine, np.stack((i, j, k), axis=1))
    return pl.DataFrame(
        {
            "i": i.astype(np.uint8),
            "j": j.astype(np.uint8),
            "k": k.astype(np.uint8),
            f_stem: values.astype(dtype),
        },
    )


def res_to_sig(f: Path) -> nb.Nifti1Image:
    tdof = nb.load(f.with_name("tdof_t1.nii.gz"))
    res4d = np.square(nb.load(f).get_fdata()).sum(axis=-1)
    sigmasquared = res4d / tdof.get_fdata()
    return nb.Nifti1Image(sigmasquared, affine=tdof.affine, header=tdof.header)


def to_polars(f: Path, dtype: np.dtype = np.float64) -> pl.DataFrame:
    if "res4d" in f.name:
        n = res_to_sig(f)
        f_stem = "sigmasquareds"
    else:
        n = nb.load(f)
        f_stem = _remove_niigz(f)
    assert len(n.shape) == 3

    return to_polars0(f_stem, n, dtype)


def join_list(
    toprocess: list[pl.DataFrame], on: list[str] = ["i", "j", "k"]
) -> pl.DataFrame:
    first = toprocess[0]
    for remaining in toprocess[1:]:
        first = first.join(remaining, on=on)
    return first


@prefect.Task
def _write(sub: Path, out: Path) -> None:
    if not out.exists():
        print(f"working on {sub=}")
        first_cope = True
        load_anats = True
        for task in sub.glob("task=*"):
            task_ = re.findall(r"task=tfMRI_(\w+)", str(task))
            assert len(task_) == 1
            for cope in task.glob("cope=*"):
                first = True
                cope_ = re.findall(r"cope=cope(\d+)", str(cope))
                assert len(cope_) == 1

                for key, dtype in FILES.items():
                    newd = to_polars(cope / key, dtype)

                    if first:
                        d = newd.clone()
                        first = False
                        target = nb.load(cope / key)

                        if load_anats:
                            anats = []
                            for anat in ANAT:
                                anatfile = sub / anat
                                resampled = processing.resample_from_to(
                                    nb.load(anatfile), target, order=0
                                )
                                anats.append(
                                    to_polars0(
                                        _remove_niigz(anatfile), resampled, np.int16
                                    )
                                )
                            anatout = join_list(anats)
                            load_anats = False

                    else:
                        d = d.join(newd, on=["i", "j", "k"])

                nrows = d.shape[0]
                d = d.with_columns(
                    pl.Series("cope", np.repeat(cope_[0], nrows).astype(np.uint8)),
                    pl.Series("task", [task_[0]] * nrows),
                )
                if first_cope:
                    dout = d.clone()
                    first_cope = False
                else:
                    dout.vstack(d, in_place=True)

        if first_cope:
            print(f"something strange while looking at {sub=}")
            return
        dout = dout.join(anatout, on=["i", "j", "k"])
        out.mkdir(parents=True)
        dout.write_parquet(out / "part-0.parquet")


@prefect.task
def get_sub(sub: Path) -> str:
    sub_ = re.findall(r"sub=(\d+)", str(sub))
    assert len(sub_) == 1
    return sub_[0]


@prefect.flow(
    # task_runner=SequentialTaskRunner
    task_runner=DaskTaskRunner(
        cluster_kwargs={
            "n_workers": 8,
            "threads_per_worker": 1,
            "dashboard_address": None,
        }
    )
)
def _main(src: Path, out: Path = Path("out")) -> None:
    for sub in src.glob("sub=*"):
        s = get_sub.submit(sub).result()
        _write.submit(sub, out=out / f"sub={s}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("src", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()

    _main(src=args.src, out=args.out)
