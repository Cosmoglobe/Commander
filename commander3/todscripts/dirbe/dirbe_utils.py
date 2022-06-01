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

DETECTORS = range(1,11)
BANDS_LABELS = ("A", "B", "C")
WAVELENGHTS = (1.25, 2.2, 3.5, 4.9, 12, 25, 60, 100, 140, 240)


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
            wavelength = WAVELENGHTS[int(detector)-1]
            band_label = f"dirbe_{int(detector):02}{band}_{wavelength}um" if band else f"dirbe_{int(detector):02}A_{wavelength}um"
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
    NSIDE = 128
    LMAX = 3 * NSIDE
    N_ALMS = LMAX**2 + 2*LMAX +1
    DEFAULT_BEAM = np.zeros((3, N_ALMS)) # Update this with actual beams

    beams: dict[str, NDArray[np.floating]] = {}
    for detector in range(1,11):
        wavelength = WAVELENGHTS[int(detector)-1]
        if detector <= 3:
            for band in BANDS_LABELS:
                beams[f"dirbe_{int(detector):02}{band}_{wavelength}um"] = DEFAULT_BEAM
        else:
            beams[f"dirbe_{int(detector):02}A_{wavelength}um"] = DEFAULT_BEAM

    return beams


def get_dirbe_sidelobes() -> dict[str, NDArray[np.floating]]:
    """
    Returns a dictionary mapping the DIRBE bands to sidelobes.
    NOTE: We dont have dirbe sidelobes so we just returns a sequence of 0's.
    """
    NSIDE = 128
    LMAX = 3 * NSIDE
    N_ALMS = LMAX**2 + 2*LMAX +1
    DEFAULT_SIDELOBES = np.zeros((3,N_ALMS)) # Update this with actual sidelobes

    sidelobes: dict[str, NDArray[np.floating]] = {}
    for detector in range(1,11):
        wavelength = WAVELENGHTS[int(detector)-1]
        if detector <= 3:
            for band in BANDS_LABELS:
                sidelobes[f"dirbe_{int(detector):02}{band}_{wavelength}um"] = DEFAULT_SIDELOBES
        else:
            sidelobes[f"dirbe_{int(detector):02}A_{wavelength}um"] = DEFAULT_SIDELOBES

    return sidelobes


if __name__ == "__main__":
    print(get_dirbe_beams().keys())