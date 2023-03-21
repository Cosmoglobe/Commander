from __future__ import annotations

from typing import TYPE_CHECKING
import sys
import numba
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

if TYPE_CHECKING:
    from ...python.commander_tools.tod_tools.commander_tod import (
        commander_tod,
    )

sys.path.insert(0, "/mn/stornext/d16/cmbco/bp/metins/Commander/commander3/python")
from commander_tools.tod_tools.commander_tod import commander_tod
import dirbe_utils

TEMP_OUTPUT_PATH = "/mn/stornext/d16/cmbco/bp/metins/dirbe/data/"
DATA_DIR = "/mn/stornext/d16/cmbco/bp/metins/dirbe/data/"
PATH_TO_HDF5_FILES = "/mn/stornext/d16/cmbco/bp/gustavbe/master/dirbe_hdf5_files/"

NSIDE = 128


@numba.njit
def accumuate_tods(emission, pixels, tods):
    for i in range(len(tods)):
        emission[pixels[i]] += tods[i]

    return emission


def main() -> None:

    tods_total = np.zeros(hp.nside2npix(NSIDE))
    pix_total = np.zeros_like(tods_total)

    comm_tod = commander_tod(TEMP_OUTPUT_PATH, "")
    N_CHUNKS = 285
    for chunk in range(284, 285):
    # for chunk in range(1, N_CHUNKS + 1):
        chunk_label = str(chunk).zfill(6)
        freq=f"DIRBE_06_25um"
        comm_tod.init_file(freq, "")
        tods = comm_tod.load_field(f'{chunk_label}/06_A/tod').astype("float")[()]
        flags = comm_tod.load_field(f'{chunk_label}/06_A/flag').astype('int')
        pix = comm_tod.load_field(f'{chunk_label}/06_A/pix').astype('int')[()]
        print(flags)
        print(tods)
        exit()
        condition = tods > 0
        print(chunk, condition.sum())
        filtered_tods = tods[condition]
        filtered_pix = pix[condition]

        tods_total = accumuate_tods(
            tods_total,
            filtered_pix,
            filtered_tods,
        )

        unique_pix, count = np.unique(filtered_pix, return_counts=True)

        pix_total[unique_pix] += count

    mask = pix_total > 0
    tods_total[mask] /= pix_total[mask]
    tods_total[~mask] = hp.UNSEEN

    hp.mollview(tods_total)
    hp.mollview(tods_total, norm="hist")
    plt.show()

if __name__ == "__main__":
    main()
    # test_smoothing()