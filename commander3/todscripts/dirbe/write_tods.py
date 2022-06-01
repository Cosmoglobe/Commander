from __future__ import annotations

from typing import TYPE_CHECKING
import sys

if TYPE_CHECKING:
    from ...python.commander_tools.tod_tools.commander_tod import (
        commander_tod,
    )

sys.path.insert(0, "/mn/stornext/u3/duncanwa/Commander/commander3/python")
from commander_tools.tod_tools.commander_tod import commander_tod
import dirbe_utils

import numpy as np
from numpy.typing import NDArray
import multiprocessing


TEMP_OUTPUT_PATH = "/mn/stornext/d16/cmbco/bp/metins/dirbe/data/"

NUM_PROCS = 1

def write_dirbe_commmander_tods(output_path: str, version: int) -> None:
    """Writes DIRBE time-ordered data to h5 files."""

    pool = multiprocessing.Pool(processes=NUM_PROCS)
    manager = multiprocessing.Manager()
    dicts = {"DIRBE": manager.dict()}

    comm_tod = commander_tod(output_path, "DIRBE", version, dicts)
    print(comm_tod)


def write_chunk(comm_tod: commander_tod, wavelength: float, chunk_idx: int) -> None:
    """Writes a single chunk of tod to file."""






def main() -> None:
    version = 1
    write_dirbe_commmander_tods(output_path=TEMP_OUTPUT_PATH, version=version)


if __name__ == "__main__":
    main()
