from __future__ import annotations
from typing import Sequence

from astropy.io import fits
import astropy.units as u
import numpy as np
from numpy.typing import NDArray
import re 


DIRBE_SKYMAP_INFO = "/mn/stornext/d16/cmbco/ola/dirbe/auxdata/DIRBE_SKYMAP_INFO.FITS"
DIRBE_BANDPASSES = "/mn/stornext/d16/cmbco/ola/dirbe/auxdata/bandpass/DIRBE_SYSTEM_SPECTRAL_RESPONSE_TABLE.ASC"
DIRBE_BEAM = "/mn/stornext/d16/cmbco/ola/dirbe/auxdata/beams/DIRBE_BEAM_CHARACTERISTICS_P3B.ASC"


def pix_to_lonlat(
    pixels: int | Sequence[int] | NDArray[np.integer],
) -> tuple[float, float] | NDArray[np.integer]:
    """Returns the lon lat pairs for DIRBE QUADCUBE pixels indices."""

    with fits.open(DIRBE_SKYMAP_INFO) as file:
        latitudes = file[1].data["ELAT-CSC"][pixels]
        longitudes = file[1].data["ELON-CSC"][pixels]

    if np.isscalar(pixels):
        return latitudes, longitudes

    return np.column_stack((latitudes, longitudes))


def normalize_dirbe_bandpasses(
    bandpasses: NDArray[np.floating],
) -> NDArray[np.floating]:
    """Normalizes the dirbe_bandpasses so that the sum of the weights is 1."""

    return bandpasses / np.expand_dims(bandpasses.sum(axis=1), axis=1)


def get_dirbe_bandpass(normalized: bool = True) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
    """Returns the DIRBE bandpasses."""

    wavelengths, *bandpasses = np.loadtxt(DIRBE_BANDPASSES, skiprows=15, unpack=True)

    wavelengths = np.asarray(wavelengths)
    bandpasses = np.asarray(bandpasses)

    if normalized:
        return wavelengths, normalize_dirbe_bandpasses(bandpasses)

    return wavelengths, bandpasses

def get_dirbe_fwhm() -> dict[str, float]:
    """Returns a dictionary mapping the DIRBE bands to FWHM in radians."""

    ROWS_TO_SKIP = 19
    FWHM_COL = 4

    fwhms: dict[str, float] = {}
    with open(DIRBE_BEAM, "r") as file:
        for line in file.readlines()[ROWS_TO_SKIP:]:
            cols = line.split()
            detector, band = re.match(r"(\d+)([A-C]?)", cols[0], re.I).groups()
            band_label = f"{int(detector):02}{band}" if band else f"{int(detector):02}A"
            if len(cols) > FWHM_COL: 
                fwhm = np.sqrt(float(cols[FWHM_COL])) * u.rad
                fwhm_arcmin = fwhm.to(u.arcmin).value
                fwhms[band_label] = np.round(fwhm_arcmin, 2)

    return fwhms


def get_dirbe_beams() -> dict[str, NDArray[np.floating]]:
    """
    Returns a dictionary mapping the DIRBE bands to beams.
    NOTE: Currently only returns a sequence of 0's.
    """

    BANDS = ("A", "B", "C")
    DEFAULT_BEAM = np.zeros((3,10000)) # Update this with actual beams

    beams: dict[str, NDArray[np.floating]] = {}
    for detector in range(1,11):
        if detector <= 3:
            for band in BANDS:
                beams[f"{detector:02}{band}"] = DEFAULT_BEAM
        else:
            beams[f"{detector:02}A"] = DEFAULT_BEAM

    return beams


def get_dirbe_sidelobes() -> dict[str, NDArray[np.floating]]:
    """
    Returns a dictionary mapping the DIRBE bands to sidelobes.
    NOTE: We dont have dirbe sidelobes so we just returns a sequence of 0's.
    """

    BANDS = ("A", "B", "C")
    DEFAULT_SIDELOBES = np.zeros((3,10000)) # Update this with actual beams

    sidelobes: dict[str, NDArray[np.floating]] = {}
    for detector in range(1,11):
        if detector <= 3:
            for band in BANDS:
                sidelobes[f"{detector:02}{band}"] = DEFAULT_SIDELOBES
        else:
            sidelobes[f"{detector:02}A"] = DEFAULT_SIDELOBES

    return sidelobes


if __name__ == "__main__":
    print(get_dirbe_sidelobes()["01A"])

