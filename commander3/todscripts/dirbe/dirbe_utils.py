from __future__ import annotations
from typing import Sequence

from astropy.io import fits
import astropy.units as u
from astropy.time import Time, TimeDelta
from datetime import datetime
import numpy as np
from numpy.typing import NDArray
import re 


DIRBE_SKYMAP_INFO = "/mn/stornext/d16/cmbco/ola/dirbe/auxdata/DIRBE_SKYMAP_INFO.FITS"
DIRBE_BANDPASSES = "/mn/stornext/d16/cmbco/ola/dirbe/auxdata/bandpass/DIRBE_SYSTEM_SPECTRAL_RESPONSE_TABLE.ASC"
DIRBE_BEAM = "/mn/stornext/d16/cmbco/ola/dirbe/auxdata/beams/DIRBE_BEAM_CHARACTERISTICS_P3B.ASC"
DIRBE_POS_PATH = "/mn/stornext/d16/cmbco/ola/dirbe/auxdata/position/"
DIRBE_POS_FILES = [
        "dmr_anc_spcl_89328_89356.txt",
        "dmr_anc_89356_90021.txt", 
        "dmr_anc_90021_90051.txt",  
        "dmr_anc_90081_90111.txt",  
        "dmr_anc_90111_90141.txt",  
        "dmr_anc_90171_90206.txt",  
        "dmr_anc_90206_90236.txt",  
        "dmr_anc_90266_90296.txt",  
        "dmr_anc_90296_90326.txt",  
        "dmr_anc_90356_91021.txt",
        "dmr_anc_90051_90081.txt", 
        "dmr_anc_90141_90171.txt",  
        "dmr_anc_90236_90266.txt",  
        "dmr_anc_90326_90356.txt",
    ]

DETECTORS = range(1,11)
DIRBE_START_DATE = Time(datetime(1989, 12, 11))
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
            band_label = f"{int(detector):02}_{band}" if band else f"{int(detector):02}_A"
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
        if detector <= 3:
            for band in BANDS_LABELS:
                beams[f"{detector:02}_{band}"] = DEFAULT_BEAM
        else:
            beams[f"{detector:02}_A"] = DEFAULT_BEAM

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
        if detector <= 3:
            for band in BANDS_LABELS:
                sidelobes[f"{detector:02}_{band}"] = DEFAULT_SIDELOBES
        else:
            sidelobes[f"{detector:02}_A"] = DEFAULT_SIDELOBES

    return sidelobes


def get_dmrfile_mjd_times(filename: str) -> NDArray[np.floating]:
    """Duncans script dont really know what goes on here."""

    data = np.loadtxt(DIRBE_POS_PATH + filename)
    adt = data[:,:2]
    i4max = 4.294967296e9
    t_adt = np.zeros_like(adt[:,0])
    ind = adt[:,0] >= 0
    t_adt[ind] = adt[ind,0] + i4max*adt[ind,1]
    ind = adt[:,0] < 0
    t_adt[ind] = i4max + adt[ind,0] + i4max*adt[ind,1]
    t_adt = t_adt*(100*u.ns)
    t_adt_d = t_adt.to('day').value # in MJD

    return t_adt_d

def get_dmrfile_positions(filename: str) -> NDArray[np.floating]:
    return np.loadtxt(DIRBE_POS_PATH + filename, usecols=[5,6,7])


def dirbe_day_to_dmr_day(days: int) -> datetime:
    """Returns the dmr_day from dirbe_day"""

    dirbe_date = DIRBE_START_DATE + TimeDelta(days - 1, format="jd")

    # datetime doesnt support leap seconds so we need to do some tricks
    dirbe_date_iso = dirbe_date.to_value(format="iso",subfmt="date_hm")

    return datetime.strptime(dirbe_date_iso, "%Y-%m-%d %H:%M")


def get_dmrfile_datetimes(filename: str) -> tuple[datetime, datetime]:
    """Returns datetime objects for the start and stop days in the dmr files."""

    regex_groups = re.search(r"(\d{2})(\d+)_(\d{2})(\d+)", filename).groups()
    year_one, yday_one, year_two, yday_two = regex_groups

    start_date = datetime.strptime(f'19{year_one} {yday_one}', '%Y %j')
    stop_date = datetime.strptime(f'19{year_two} {yday_two}', '%Y %j')
    
    return start_date, stop_date


def test_naming() -> None:
    import h5py
    with h5py.File("test.h5", "w") as file:
        for detector in DETECTORS:
            file.create_group(f"{detector:02}_{WAVELENGHTS[detector - 1]}um")
            for band in BANDS_LABELS:
                file.create_group(f"{detector:02}_{band}")
                if detector > 3:
                    break

if __name__ == "__main__":
    ...