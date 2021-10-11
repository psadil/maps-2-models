import os
from typing import Union, Optional

import nimare
from nimare.meta import ibma

def do_imba(
  dset_fname: Union[str, bytes, os.PathLike],
  in_dir: Union[str, bytes, os.PathLike],
  out_dir: Union[str, bytes, os.PathLike] = os.getcwd(),
  prefix: str = "prefix"
  ) -> Union[str, bytes, os.PathLike]:
          
    dset = nimare.dataset.Dataset.load(dset_fname)
    dset.update_path(in_dir)
    meta = ibma.Stouffers(use_sample_size=True, resample=False)
    meta.fit(dset)
    meta.results.save_maps(output_dir=out_dir, prefix=prefix)
    meta.save(os.path.join(out_dir, f"{prefix}_dset.pklz"))

    return os.path.join(out_dir, f"{prefix}_z.nii.gz")
