from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass

from astropy.io import fits
import numpy as np

CIO_PATH = Path("/mn/stornext/d16/cmbco/ola/dirbe/cio")
MIN_OBS_DAY = 89345
MAX_OBS_DAY = 90264

@dataclass
class ChunkData:
    tods: np.ndarray
    pix: np.ndarray
    time: np.ndarray

    def from_fits(self, file: Path) -> ChunkData:
        if not file.exists():
            raise FileNotFoundError(f"File not found: {file}")

        with fits.open(file) as hdul:
            data = hdul[1].data
            pix9 = data["Pixel_no"].astype(np.int64)
            pix9_sub = data["PSubPos"].astype(np.int64)
            pix9_sub_sub = data["PSbSbPos"].astype(np.int64)
            pix15 = 4096 * pix9 + 16*pix9_sub + pix9_sub_sub


def get_res_15_pixels(yday: int) -> np.ndarray:
    file = CIO_PATH / f"DIRBE_CIO_P3B_{yday:5}.FITS"
    if not file.exists():
        raise FileNotFoundError(f"File not found: {file}")

    with fits.open(file) as hdul:
        data = hdul[1].data
        print(data[""])
        pix9 = data["Pixel_no"].astype(np.int64)
        pix9_sub = data["PSubPos"].astype(np.int64)
        pix9_sub_sub = data["PSbSbPos"].astype(np.int64)
        pix15 = 4096 * pix9 + 16*pix9_sub + pix9_sub_sub

    return pix15

if __name__ == "__main__":
    pix = get_res_15_pixels(89345)
    print(pix)


def get_time_ordered_inds(yday: int) -> np.ndarray:
    pix = get_res_15_pixels(yday)
    inds = np.argsort(pix)
    return inds