import os
from typing import Union, Optional

import pandas as pd

import nimare
from nimare.meta.cbma import ALE
from nimare.transforms import ImagesToCoordinates

def as_dict(d):
  ds = {
    d.exp: {
      'contrasts': {
        "1": {
          "images": {
            'z': d.z
            },
            "metadata": {
              "sample_sizes": [d.n_sub]
              } 
            }
          }
        }
      }

  return ds

def do_ale(
  d: pd.DataFrame,
  in_dir: Union[str, bytes, os.PathLike],
  out_dir: Union[str, bytes, os.PathLike] = os.getcwd(),
  prefix: Optional[str] = None,
  mask: Union[str, bytes, os.PathLike] = os.path.join(os.getenv("FSLDIR"), "data", "standard", "MNI152_T1_2mm_brain_mask.nii.gz")
) -> Union[str, bytes, os.PathLike]:

    x = {}
    for row in d.itertuples():
      x.update(as_dict(row))
          
    dset = nimare.dataset.Dataset(x, mask=mask)
    dset.update_path(new_path=in_dir)
    coord_replace = ImagesToCoordinates(merge_strategy="replace", z_threshold=1, remove_subpeaks=True)
    dset = coord_replace.transform(dset)   
    dset.save(os.path.join(out_dir, f"{prefix}_dset.pklz"))
    ale = ALE()

    cbma = ale.fit(dset)
    cbma.save_maps(output_dir=out_dir, prefix=prefix)

    return os.path.join(out_dir, f"{prefix}_z.nii.gz")

