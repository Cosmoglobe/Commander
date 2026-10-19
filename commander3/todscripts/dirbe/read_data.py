from __future__ import annotations
from pathlib import Path
from functools import cache
from dataclasses import dataclass
import quadcube

from astropy.io import fits
import astropy.units as u
import numpy as np

CIO_PATH = Path("/mn/stornext/d16/cmbco/ola/dirbe/cio")
BEAM_FILE = Path("/mn/stornext/d16/cmbco/ola/dirbe/DIRBE_BEAM_CHARACTERISTICS_P3B.ASC")
DETECTORS = [
    "1A",
    "1B",
    "1C",
    "2A",
    "2B",
    "2C",
    "3A",
    "3B",
    "3C",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
]

ROWS_IN_BEAM_FILE_TO_SKIP = 19
MIN_OBS_DAY = 89345
MAX_OBS_DAY = 90264

YDAYS = np.concatenate([np.arange(89345, 89366), np.arange(90001, 90265)])
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
            pix15 = 4096 * pix9 + 16 * pix9_sub + pix9_sub_sub


def get_cio_data(yday: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    file = CIO_PATH / f"DIRBE_CIO_P3B_{yday:5}.FITS"
    if not file.exists():
        raise FileNotFoundError(f"File not found: {file}")

    with fits.open(file) as hdul:
        data = hdul[1].data

        pix9 = data["Pixel_no"].astype(np.int64)
        pix9_sub = data["PSubPos"].astype(np.int64)
        pix9_sub_sub = data["PSbSbPos"].astype(np.int64)
        pix15 = 4096 * pix9 + 16 * pix9_sub + pix9_sub_sub

        attack_vecs = data["AttackV"].astype(np.float64)
        natv = attack_vecs / np.expand_dims(np.linalg.norm(attack_vecs, axis=1), axis=1)

        time = data["Time"].astype(np.float64)

    return pix15, natv, time


@dataclass
class DetectorBeamData:
    solid_angle: float
    solid_angle_error: float
    beam: dict[range, tuple[float, float]]


@cache
def get_beam_data() -> dict[str, DetectorBeamData]:
    data: dict[str, DetectorBeamData] = {}
    visited = set()
    current_det_data = {}
    with open(BEAM_FILE, "r") as f:
        lines = iter(f.readlines()[ROWS_IN_BEAM_FILE_TO_SKIP:])
        while True:
            try:
                line = next(lines)
            except StopIteration:
                data[split[0]] = DetectorBeamData(**current_det_data)
                break
            split = line.split()
            current_det = split[0]
            if current_det not in visited:
                visited.add(current_det)
                current_det_data = {
                    "solid_angle": float(split[4]),
                    "solid_angle_error": float(split[5]),
                    "beam": {},
                }
                if visited != set():
                    data[split[0]] = DetectorBeamData(**current_det_data)
            if current_det in visited:
                start, stop = split[1].split("-")
                current_det_data["beam"][range(int(start), int(stop))] = (
                    u.Quantity(float(split[2]), u.arcmin).to_value(u.arcmin),
                    u.Quantity(float(split[3]), u.arcmin).to_value(u.arcmin),
                )

    return data


def get_unit_vectors():
    for yday in range(MIN_OBS_DAY, MAX_OBS_DAY + 1):
        pix, natv, time = get_cio_data(yday)
        los = quadcube.pix2vec(pix)


        beam_data = get_beam_data()
        for band in DETECTORS:
            for r in beam_data[band].beam:
                if yday in r:
                    isco, xsco = beam_data[band].beam[r]
                    break
            else: 
                raise ValueError(f"No beam data for {yday} and {band}")

            new_los = np.zeros_like(los)
            new_los[0] = los[0] + isco * natv[:, 0] + xsco * (los[1] * natv[:, 2] - los[2] * natv[:, 1])
            new_los[1] = los[1] + isco * natv[:, 1] + xsco * (los[2] * natv[:, 0] - los[0] * natv[:, 2])
            new_los[2] = los[2] + isco * natv[:, 2] + xsco * (los[0] * natv[:, 1] - los[1] * natv[:, 0])

            import matplotlib.pyplot as plt
            inds = np.argsort(time)
            sorted_now_los = new_los[:, inds]
            plt.plot(sorted_now_los[0, :10000])
            plt.show()
            exit()

if __name__ == "__main__":
    get_unit_vectors()
