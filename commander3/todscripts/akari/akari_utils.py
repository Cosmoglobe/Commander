from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence

from functools import cache

from astropy.io import fits
import astropy.units as u
from astropy.time import Time, TimeDelta
from datetime import datetime
import numpy as np
from numpy.typing import NDArray
import re
from scipy.interpolate import interp1d
import healpy as hp
from astropy.coordinates import (
    solar_system_ephemeris,
    get_body,
    HeliocentricMeanEcliptic,
)

DIRBE_SKYMAP_INFO = "/mn/stornext/d16/cmbco/ola/dirbe/auxdata/DIRBE_SKYMAP_INFO.FITS"
DIRBE_BANDPASSES = "/mn/stornext/d16/cmbco/ola/dirbe/auxdata/bandpass/DIRBE_SYSTEM_SPECTRAL_RESPONSE_TABLE.ASC"
DIRBE_BEAM = (
    "/mn/stornext/d16/cmbco/ola/dirbe/auxdata/beams/DIRBE_BEAM_CHARACTERISTICS_P3B.ASC"
)
DIRBE_POS_PATH = "/mn/stornext/d16/cmbco/ola/dirbe/auxdata/position/"
BEAM_FILE = "/mn/stornext/d16/cmbco/ola/dirbe/DIRBE_BEAM_CHARACTERISTICS_P3B.ASC"
BANDPASS_PATH = "/mn/stornext/d5/data/metins/dirbe/data/"

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

BANDS = range(1, 5)
DIRBE_START_DATE = Time(datetime(2006, 4, 1))   # April 1, 2006 Date uncertain
DETECTOR_LABELS = ("A", "B", "C")
WAVELENGTHS = (65, 90, 140, 160)
DETECTORS = [
        'N60',         # Number of pixels (40) 25.28 Hz fsamp
        'WIDE-S',      #                  (60) 25.28 Hz fsamp
        'WIDE-L',      #                  (45)  16.86 fsamp
        'N160'         #                  (30)  16.86 fsamp
]

BAND_TO_WAVELENGTH: dict[int, float] = {
    1: 1.25,
    2: 2.2,
    3: 3.5,
    4: 4.9,
    5: 12,
    6: 25,
    7: 60,
    8: 100,
    9: 140,
    10: 240,
}

ROWS_IN_BEAM_FILE_TO_SKIP = 19

PLANET_RADII = {
    "moon": 10,
    "mercury": 1,
    "venus": 1,
    "mars": 1,
    "jupiter": 1.5,
    "saturn": 1,
    "uranus": 1,
    "neptune": 1,
}

SIGMA_0 = {
        1: 40
}


def get_planet_interps(time_delta: TimeDelta) -> dict[str, dict[str, interp1d]]:
    times = np.arange(datetime(1989, 6, 1), datetime(1991, 1, 1), time_delta).astype(
        datetime
    )
    astropy_times = Time(times, format="datetime", scale="utc")
    interpolaters = {}
    rotator = hp.Rotator(coord=["E", "G"])
    with solar_system_ephemeris.set("de432s"):
        for body_name in PLANET_RADII:
            interpolaters[body_name] = {}
            body = get_body(body_name, astropy_times).transform_to(
                "geocentricmeanecliptic"
            )
            lon, lat = rotator(body.lon.value, body.lat.value, lonlat=True)
            interpolaters[body_name]["lon"] = interp1d(astropy_times.mjd, lon)
            interpolaters[body_name]["lat"] = interp1d(astropy_times.mjd, lat)

    return interpolaters


@cache
def band_to_bit(band: int) -> int:
    """Returns the bit corresponding to the band."""
    return {
        1: 0,
        2: 3,
        3: 6,
        4: 9,
        5: 10,
        6: 11,
        7: 12,
        8: 13,
        9: 14,
        10: 15,
    }[band]


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


@cache
def get_akari_fwhm() -> dict[str, float]:
    """Returns a dictionary mapping the Akari bands to FWHM in radians."""

    ROWS_TO_SKIP = 19
    FWHM_COL = 4

    fwhms: dict[str, float] = {}
    with open(DIRBE_BEAM, "r") as file:
        for line in file.readlines()[ROWS_TO_SKIP:]:
            cols = line.split()
            detector, band = re.match(r"(\d+)([A-C]?)", cols[0], re.I).groups()
            band_label = (
                f"{int(detector):02}_{band}" if band else f"{int(detector):02}_A"
            )
            if len(cols) > FWHM_COL:
                fwhm = np.sqrt(float(cols[FWHM_COL])) * u.rad
                fwhm_arcmin = fwhm.to(u.arcmin).value
                fwhms[band_label] = np.round(fwhm_arcmin, 2)

    return fwhms


@cache
def get_akari_beams() -> dict[str, NDArray[np.floating]]:
    """
    Returns a dictionary mapping the DIRBE bands to beams.
    NOTE: Currently only returns a sequence of 0's.
    """
    NSIDE = 128
    LMAX = 3 * NSIDE
    N_ALMS = LMAX**2 + 2 * LMAX + 1
    DEFAULT_BEAM = np.zeros((3, N_ALMS))  # Update this with actual beams

    beams: dict[str, NDArray[np.floating]] = {}
    for detector in range(1, 11):
        if detector <= 3:
            for band in DETECTOR_LABELS:
                beams[f"{detector:02}_{band}"] = DEFAULT_BEAM
        else:
            beams[f"{detector:02}_A"] = DEFAULT_BEAM

    return beams


@cache
def get_akari_sidelobes() -> dict[str, NDArray[np.floating]]:
    """
    Returns a dictionary mapping the DIRBE bands to sidelobes.
    NOTE: We dont have dirbe sidelobes so we just returns a sequence of 0's.
    """
    NSIDE = 128
    LMAX = 3 * NSIDE
    N_ALMS = LMAX**2 + 2 * LMAX + 1
    DEFAULT_SIDELOBES = np.zeros((3, N_ALMS))  # Update this with actual sidelobes

    sidelobes: dict[str, NDArray[np.floating]] = {}
    for detector in range(1, 11):
        if detector <= 3:
            for band in DETECTOR_LABELS:
                sidelobes[f"{detector:02}_{band}"] = DEFAULT_SIDELOBES
        else:
            sidelobes[f"{detector:02}_A"] = DEFAULT_SIDELOBES

    return sidelobes


@cache
def get_dmrfile_mjd_times(filename: str) -> NDArray[np.floating]:
    """Duncans script dont really know what goes on here."""

    data = np.loadtxt(DIRBE_POS_PATH + filename)
    adt = data[:, :2]
    i4max = 4.294967296e9
    t_adt = np.zeros_like(adt[:, 0])
    ind = adt[:, 0] >= 0
    t_adt[ind] = adt[ind, 0] + i4max * adt[ind, 1]
    ind = adt[:, 0] < 0
    t_adt[ind] = i4max + adt[ind, 0] + i4max * adt[ind, 1]
    t_adt = t_adt * (100 * u.ns)
    t_adt_d = t_adt.to("day").value  # in MJD

    return t_adt_d


@cache
def get_dmrfile_positions(filename: str) -> NDArray[np.floating]:
    return np.loadtxt(DIRBE_POS_PATH + filename, usecols=[5, 6, 7])


@cache
def dirbe_day_to_dmr_day(days: int) -> datetime:
    """Returns the dmr_day from dirbe_day"""

    dirbe_date = DIRBE_START_DATE + TimeDelta(days - 1, format="jd")

    # datetime doesnt support leap seconds so we need to do some tricks
    dirbe_date_iso = dirbe_date.to_value(format="iso", subfmt="date_hm")

    return datetime.strptime(dirbe_date_iso, "%Y-%m-%d %H:%M")


@cache
def get_dmrfile_datetimes(filename: str) -> tuple[datetime, datetime]:
    """Returns datetime objects for the start and stop days in the dmr files."""

    regex_groups = re.search(r"(\d{2})(\d+)_(\d{2})(\d+)", filename).groups()
    year_one, yday_one, year_two, yday_two = regex_groups

    start_date = datetime.strptime(f"19{year_one} {yday_one}", "%Y %j")
    stop_date = datetime.strptime(f"19{year_two} {yday_two}", "%Y %j")

    return start_date, stop_date


def test_naming() -> None:
    import h5py

    with h5py.File("test.h5", "w") as file:
        for detector in BANDS:
            file.create_group(f"{detector:02}_{WAVELENGHTS[detector - 1]}um")
            for band in DETECTOR_LABELS:
                file.create_group(f"{detector:02}_{band}")
                if detector > 3:
                    break


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
                if current_det in {"9", "10"}:
                    start, stop = "89345", "90260"
                # Hack to invlude final 5 days fo data...
                if stop == "90260":
                    stop = "90266"
                current_det_data["beam"][range(int(start), int(stop))] = (
                    u.Quantity(float(split[2]), u.arcmin).to_value(u.rad),
                    u.Quantity(float(split[3]), u.arcmin).to_value(u.rad),
                )

    return data


def get_sat_and_earth_pos(yday: int, time: float) -> tuple[np.ndarray, np.ndarray]:
    """dmr_cio_91206-91236.fits contains data from day 206 of 1991 to day 236 of 1991)."""

    ranges_times = [re.findall(r"\d+", file) for file in DIRBE_POS_FILES]
    ranges = [range(*[int(time) for time in times]) for times in ranges_times]
    for idx, r in enumerate(ranges):
        if yday in r:
            dmr_file = DIRBE_POS_FILES[idx]
            break
    else:
        raise FileNotFoundError(f"could not find file containing sat pos for {yday=}")

    dmr_times = get_dmrfile_mjd_times(dmr_file)
    dmr_positions = get_dmrfile_positions(dmr_file)

    interpolator_sat_x = interp1d(
        dmr_times,
        dmr_positions[:, 0],
        fill_value="extrapolate",
    )
    interpolator_sat_y = interp1d(
        dmr_times,
        dmr_positions[:, 1],
        fill_value="extrapolate",
    )
    interpolator_sat_z = interp1d(
        dmr_times,
        dmr_positions[:, 2],
        fill_value="extrapolate",
    )

    pos_x = interpolator_sat_x(time)
    pos_y = interpolator_sat_y(time)
    pos_z = interpolator_sat_z(time)

    celestial_sat_pos = u.Quantity([pos_x, pos_y, pos_z], u.m).to(u.AU).value

    rotator = hp.Rotator(coord=["C", "E"])
    geocentric_ecl_sat_pos = u.Quantity(rotator(celestial_sat_pos), u.AU).transpose()

    earth_pos = get_body("earth", Time(time, format="mjd")).transform_to(
        HeliocentricMeanEcliptic
    )
    earth_pos = earth_pos.cartesian.xyz.to(u.AU).transpose()

    ecl_sat_pos = earth_pos + geocentric_ecl_sat_pos

    return ecl_sat_pos.value, earth_pos.value


@cache
def get_const_scalars(band: int) -> NDArray[np.floating]:
    """Used in V14 -> and out"""
    # SIGMA0 = u.Quantity(
    #     [2.4, 1.6, 0.9, 0.8, 0.9, 0.9, 0.9, 0.5, 32.8, 10.7], "nW/(m^2 sr)"
    # )
    # SIGMA0 *= 20
    # SIGMA0 /= u.Quantity(
    #     [59.5, 22.4, 22.0, 8.19, 13.3, 4.13, 2.32, 0.974, 0.605, 0.495], "THz"
    # )
    # SIGMA0 = SIGMA0.to_value("MJy/sr")

    TEMP_GAIN = 1
    TEMP_ALPHA = -1
    fknee = 1 / (10 * 60)

    return np.array([TEMP_GAIN, SIGMA_0[band], fknee, TEMP_ALPHA]).flatten()

@cache
def get_bandpass(band: int) -> tuple[u.Quantity[u.micron], NDArray[np.float64]]:
    bandpass_file = BANDPASS_PATH + f"/DIRBE_{band:02}_bandpass.dat"
    bandpass = np.loadtxt(bandpass_file, unpack=True)

    non_zero = np.nonzero(bandpass[1])
    bandpass = bandpass[:, non_zero[0]]
    freqs, weights = bandpass
    freqs *= u.micron

    return freqs, weights


@cache
def get_iras_factor(band: int) -> float:
    freqs, weights = get_bandpass(band)
    freq_ref = BAND_TO_WAVELENGTH[band]
    freqs = freqs.to(u.Hz, u.spectral())
    freqs = np.flip(freqs)
    weights = np.flip(weights)
    weights /= np.trapz(weights, freqs.value)
    return np.trapz(
        weights * ((freq_ref * u.micron).to(u.Hz, u.spectral()).value / freqs.value),
        freqs.value,
    )


if __name__ == "__main__":
    print([get_iras_factor(i) for i in range(1,11)])

