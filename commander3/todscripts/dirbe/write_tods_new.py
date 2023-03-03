"""
version 2: add smoothing of pixels to remove pixelization effects
version 3: rotate tods to galactic coordinates
version 4: add custom planet flags
version 5: larger flagging radii
version 6: per pointing time for accurate planet flags
version 7: redo planet flags
version 8: redo planet flags better
version 9: fix scalars, try smaller flagging radius
version 10: gapfilling so that all tods are equispaced
version 11: change smoothing to using gaussian filter and remove iras convention color correction from the tods
"""

from __future__ import annotations
from multiprocessing.managers import DictProxy
import os
import random

from typing import TYPE_CHECKING
import sys
from scipy import interpolate
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from pathlib import Path

if TYPE_CHECKING:
    from ...python.commander_tools.tod_tools.commander_tod import (
        commander_tod,
    )

sys.path.insert(0, "/mn/stornext/d5/data/metins/Commander/commander3/python")
from commander_tools.tod_tools.commander_tod import commander_tod
import dirbe_utils

from astropy.coordinates import (
    HeliocentricMeanEcliptic,
    get_body,
    solar_system_ephemeris,
)
from scipy.interpolate import interp1d
import healpy as hp
import numpy as np
from numpy.typing import NDArray
import multiprocessing
import h5py
import astropy.units as u
from astropy.time import Time

TEMP_OUTPUT_PATH = "/mn/stornext/d5/data/metins/dirbe/data/"
PATH_TO_HDF5_FILES = "/mn/stornext/d16/cmbco/bp/gustavbe/master/dirbe_hdf5_files/"
BANDPASS_PATH = Path("/mn/stornext/d5/data/metins/dirbe/data/")

NUM_PROCS = 10
NUM_CHUNKS_PER_BAND = 285
DEFAULT_NSIDE = 128
NPSI = 2048
FSAMP = 8
BAD_DATA_SENTINEL = -16375
N_CONV_BOX_POINTS = 30
NSIDE = 128

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

BODIES = {
    "moon": 1,
    "mercury": 1,
    "venus": 1,
    "mars": 1,
    "jupiter": 1,
    "saturn": 1,
    "uranus": 1,
    "neptune": 1,
}


def get_flag_sum(*flag_bits: int) -> int:
    return sum([2**bit for bit in flag_bits])


flag_bit_sum = get_flag_sum(
    0,  # radiation zone
    3,  # oa
    13,  # bad data (sentinel values)
    14,  # moon
    15,  # mercury
    16,  # venus
    17,  # mars
    18,  # jupiter
    19,  # saturn
    20,  # uranus
    21,  # neptune
)


def write_dirbe_commmander_tods(
    output_path: str,
    version: int,
    smooth_pixels: bool = False,
    nside_in: int = NSIDE,
    nside_out: int = NSIDE,
) -> None:
    """Writes DIRBE time-ordered data to h5 files."""

    random.seed()
    os.environ["OMP_NUM_THREADS"] = "1"
    pool = multiprocessing.Pool(processes=NUM_PROCS)
    manager = multiprocessing.Manager()

    times = np.arange(
        datetime(1989, 6, 1), datetime(1991, 1, 1), timedelta(hours=1)
    ).astype(datetime)
    astropy_times = Time(times, format="datetime", scale="utc")
    planet_interps = get_planet_interps(astropy_times)

    multiprocessor_manager_dicts: dict[str, DictProxy] = {}
    for idx, wavelength in enumerate(dirbe_utils.WAVELENGHTS, start=1):
        name = f"DIRBE_{idx:02}_{wavelength}um"
        if smooth_pixels:
            name += "_smoothed"
        if nside_out != nside_in:
            name += f"_nside{nside_out}"
        name += f"_V{version:02}"
        multiprocessor_manager_dicts[name] = manager.dict()

    filenames = list(multiprocessor_manager_dicts.keys())
    comm_tod = commander_tod(output_path, "", version, multiprocessor_manager_dicts)
    x = [
        [
            pool.apply_async(
                write_band,
                args=[
                    comm_tod,
                    band,
                    nside_in,
                    nside_out,
                    smooth_pixels,
                    filenames[band - 1],
                    planet_interps,
                ],
            )
            for band in dirbe_utils.BANDS
        ]
    ]

    for res1 in x:
        for res in res1:
            res.get()

    pool.close()
    pool.join()

    comm_tod.make_filelists()


def write_band(
    comm_tod: commander_tod,
    band: int,
    nside_in: int,
    nside_out: int,
    smooth: bool,
    filename: str,
    planet_interps: dict[str, dict[str, interp1d]],
) -> None:
    """Writes a single chunk of tod to file."""
    hdf5_filename = PATH_TO_HDF5_FILES + f"Phot{band:02}_{nside_in}.hdf5"
    COMMON_GROUP = "/common"
    HUFFMAN_COMPRESSION = ["huffman", {"dictNum": 1}]

    comm_tod.init_file(freq=filename, od="", mode="w")

    comm_tod.add_field(COMMON_GROUP + "/fsamp", FSAMP)

    # nside
    comm_tod.add_field(COMMON_GROUP + "/nside", [nside_out])

    # make detector names lookup
    if band <= 3:
        detector_names = ",".join(
            [f"{band:02}_{det_label}" for det_label in dirbe_utils.DETECTOR_LABELS]
        )
    else:
        detector_names = f"{band:02}_{dirbe_utils.DETECTOR_LABELS[0]},"
    comm_tod.add_field(COMMON_GROUP + "/det", np.string_(detector_names))

    common_polang = get_band_polang(band)
    comm_tod.add_field(COMMON_GROUP + "/polang", common_polang)
    comm_tod.add_attribute(COMMON_GROUP + "/polang", "index", detector_names)

    common_mbang = get_mbang(band)
    comm_tod.add_field(COMMON_GROUP + "/mbang", common_mbang)
    comm_tod.add_attribute(COMMON_GROUP + "/mbang", "index", detector_names)

    rotator = hp.Rotator(coord=["E", "G"])
    pid = 0
    iras_factor = get_iras_factor(band)
    for chunk in range(1, NUM_CHUNKS_PER_BAND + 1):
        day = chunk
        print(f"processing chunk: {chunk:03}...")
        chunk_label = f"{chunk:06}"

        mjd_times = get_chunk_mjd_times(chunk_label, hdf5_filename)
        if mjd_times.size == 0:
            continue
        fsamp_in_days = 1 / FSAMP / 86400
        inds = np.nonzero(np.diff(mjd_times) > (2 * fsamp_in_days))[0] + 1
        tods_det = []
        pixels_det = []
        flags_det = []
        psi_det = []
        scalars_det = []
        out_ang_det = []

        for det in dirbe_utils.DETECTOR_LABELS:
            tods = get_chunk_band_tods(chunk_label, hdf5_filename, det)
            tods *= iras_factor
            pixels = get_chunk_band_pixels(chunk_label, hdf5_filename, det)
            gustav_flags = get_gustav_flags(chunk_label, hdf5_filename, det)
            gustav_time = get_gustav_time(chunk_label, hdf5_filename, det)
            flags = get_chunk_band_flags(
                nside_out,
                planet_interps,
                pixels,
                tods,
                gustav_flags,
                gustav_time,
            )
            psi = get_chunk_band_psi(chunk_label, hdf5_filename, det)
            out_ang = get_out_ang()
            scalars = get_scalars(tods, flags)

            tods_det.append(tods)
            pixels_det.append(pixels)
            flags_det.append(flags)
            psi_det.append(psi)
            scalars_det.append(scalars)
            out_ang_det.append(out_ang)
            if band > 3:
                break

        start_idx = 0
        for idx in range(len(inds)):
            pid += 1
            pid_label = f"{pid:06}"
            pid_common_group = pid_label + "/common"
            stop_idx = inds[idx]

            if (stop_idx - start_idx) == 0:
                continue
            mjd_times_pid = mjd_times[start_idx:stop_idx]
            comm_tod.add_field(
                pid_common_group + "/time", [mjd_times_pid[0], 0, 0]
            )  # TODO: get SCET and OBT times
            comm_tod.add_attribute(
                pid_common_group + "/time", "index", "MJD, OBT, SCET"
            )
            ntod = len(tods_det[0][start_idx:stop_idx])
            comm_tod.add_field(pid_common_group + "/ntod", [ntod])

            sat_pos, earth_pos = get_sat_and_earth_pos(day, mjd_times_pid)
            comm_tod.add_field(pid_common_group + "/satpos", sat_pos[0])

            # Add earth position. Required for zodiacal light calculation.
            comm_tod.add_field(pid_common_group + "/earthpos", earth_pos[0])

            # add metadata
            comm_tod.add_attribute(pid_common_group + "/satpos", "index", "X, Y, Z")
            comm_tod.add_attribute(
                pid_common_group + "/satpos", "coords", "heliocentric-ecliptic"
            )
            for det_idx, det in enumerate(dirbe_utils.DETECTOR_LABELS):
                pid_data_group = f"{pid_label}/{band:02}_{det}"
                flags_pid = flags_det[det_idx][start_idx:stop_idx]
                tods_pid = tods_det[det_idx][start_idx:stop_idx]
                pixels_pid = pixels_det[det_idx][start_idx:stop_idx]
                psi_pid = psi_det[det_idx][start_idx:stop_idx]
                scalars_pid = scalars_det[det_idx]
                out_ang_pid = out_ang_det[det_idx]

                pixels_pid = smooth_and_udgrade_and_rotate_pixels(
                    pixels_pid, nside_in, nside_out, rotator, smooth
                )

                if flags_pid.size > 0:
                    comm_tod.add_field(
                        pid_data_group + "/flag", flags_pid, HUFFMAN_COMPRESSION
                    )

                if tods_pid.size > 0:
                    comm_tod.add_field(pid_data_group + "/tod", tods_pid)

                if pixels_pid.size > 0:
                    comm_tod.add_field(
                        pid_data_group + "/pix", pixels_pid, HUFFMAN_COMPRESSION
                    )

                psi_digitize_compression = [
                    "digitize",
                    {"min": 0, "max": 2 * np.pi, "nbins": NPSI},
                ]
                if psi_pid.size > 0:
                    comm_tod.add_field(
                        pid_data_group + "/psi",
                        psi_pid,
                        [psi_digitize_compression, HUFFMAN_COMPRESSION],
                    )

                comm_tod.add_field(pid_data_group + "/outP", out_ang_pid)

                comm_tod.add_field(pid_data_group + "/scalars", scalars_pid)
                comm_tod.add_attribute(
                    pid_data_group + "/scalars", "index", "gain, sigma0, fknee, alpha"
                )

                if band > 3:
                    break

            comm_tod.finalize_chunk(pid)
            start_idx = stop_idx

    comm_tod.finalize_file()


def smooth_and_udgrade_and_rotate_pixels(
    pixels: NDArray[np.int64], nside_in: int, nside_out: int, rotator: hp.Rotator, smooth: bool
) -> NDArray[np.int64]:
    """Smooth pixels and then udgrade them to the desired nside."""
    unit_vectors = np.asarray(hp.pix2vec(nside_in, pixels))
    if smooth: 
        unit_vectors = np.array(
            [
                gaussian_filter1d(unit_vector, sigma=5, mode="nearest")
                for unit_vector in unit_vectors
            ]
        )
    unit_vectors_gal = rotator(unit_vectors)
    return hp.vec2pix(nside_out, *unit_vectors_gal)


def get_band_polang(band: int) -> list[float]:
    """TODO: get correct angles"""

    if band > 3:
        return [0]
    return [0, 0, 0]


def get_mbang(band: int) -> list[float]:
    """TODO: get correct angles"""

    if band > 3:
        return [0]
    return [0, 0, 0]


def get_chunk_mjd_times(chunk_label: str, hdf5_filename: str) -> NDArray[np.floating]:
    """Gets dirbe mjd times from gustavs files."""

    with h5py.File(hdf5_filename, "r") as file:
        return file[f"{chunk_label}/A/time"][()]


def get_planet_interps(astropy_times: Time) -> dict[str, dict[str, interp1d]]:
    print("precomputing planet interpolators...")
    interpolaters = {}
    rotator = hp.Rotator(coord=["E", "G"])
    with solar_system_ephemeris.set("de432s"):
        for body_name in BODIES:
            interpolaters[body_name] = {}
            body = get_body(body_name, astropy_times).transform_to(
                "geocentricmeanecliptic"
            )
            lon, lat = rotator(body.lon.value, body.lat.value, lonlat=True)
            interpolaters[body_name]["lon"] = interp1d(astropy_times.mjd, lon)
            interpolaters[body_name]["lat"] = interp1d(astropy_times.mjd, lat)

    print("done")
    return interpolaters


def get_gustav_flags(
    chunk_label: str, hdf5_filename: str, detector: str
) -> NDArray[np.integer]:
    with h5py.File(hdf5_filename, "r") as file:
        return file[f"{chunk_label}/{detector}/flag"][()]


def get_gustav_time(
    chunk_label: str, hdf5_filename: str, detector: str
) -> NDArray[np.integer]:
    with h5py.File(hdf5_filename, "r") as file:
        return file[f"{chunk_label}/{detector}/time"][()]


def get_chunk_band_flags(
    nside_out: int,
    planet_interpolaters: dict[str, dict[str, interp1d]],
    pixels: NDArray[np.integer],
    tods: NDArray[np.floating],
    flags: NDArray[np.integer],
    time: NDArray[np.floating],
) -> NDArray[np.integer]:
    """Reprocess gustavs flags.

    Gustav added the following flags:
    - 0: radiation zone
    - 3: oa
    - 10: excess noise
    - 11: moon
    - 12: jupiter

    We remove 11 and 12 (because we add redo these flags)

    We add flags:
    - 13: bad data sentinel
    - 14-x: see `BODIES` object.

    """
    flags[tods <= BAD_DATA_SENTINEL] += int(2**13)

    for bit_idx, (body, radius) in enumerate(BODIES.items(), start=14):
        lon, lat = hp.pix2ang(nside_out, pixels, lonlat=True)
        planet_lon = planet_interpolaters[body]["lon"](time)
        planet_lat = planet_interpolaters[body]["lat"](time)
        ang_dist = hp.rotator.angdist(
            np.asarray([lon, lat]),
            np.asarray([planet_lon, planet_lat]),
            lonlat=True,
        )

        planet_indices = ang_dist <= np.deg2rad(radius)
        flags[planet_indices] += int(2 ** (bit_idx))

    return flags


def get_chunk_band_tods(
    chunk_label: str, hdf5_filename: str, detector: str
) -> NDArray[np.floating]:
    """Gets dirbe tods from gustavs files."""

    with h5py.File(hdf5_filename, "r") as file:
        return file[f"{chunk_label}/{detector}/tod"][()]


def get_ntods(chunk_label: str, hdf5_filename: str) -> int:
    """Gets dirbe tods from gustavs files."""

    with h5py.File(hdf5_filename, "r") as file:
        return len(file[f"{chunk_label}/A/tod"])


def get_chunk_band_psi(
    chunk_label: str, hdf5_filename: str, detector: str
) -> NDArray[np.floating]:
    """Gets dirbe polang from gustavs files. TODO: rotate to galactic."""

    with h5py.File(hdf5_filename, "r") as file:
        return file[f"{chunk_label}/{detector}/polang"][()]


def get_chunk_band_pixels(
    chunk_label: str, hdf5_filename: str, detector: str
) -> NDArray[np.integer]:
    """Gets dirbe pixels from gustavs files."""

    with h5py.File(hdf5_filename, "r") as file:
        pixels_ecl = file[f"{chunk_label}/{detector}/pix"][()]

    return pixels_ecl


def get_out_ang() -> NDArray[np.floating]:
    """TODO: get correct out ang maybe..?"""

    return np.zeros((2, 1))


def get_scalars(tods: NDArray[np.floating], flags) -> NDArray[np.floating]:
    """TODO: get correct gain, sigma0, fknee and alpha."""

    TEMP_GAIN = 1
    TEMP_ALPHA = -1
    fknee = 1 / (10 * 60)

    condition = np.bitwise_and(flags, flag_bit_sum) == 0
    filtered_tods = tods[condition]
    sigma0 = np.diff(filtered_tods).std() / (2**0.5)

    return np.array([TEMP_GAIN, sigma0, fknee, TEMP_ALPHA]).flatten()


def get_sat_and_earth_pos(
    day: int, dirbe_times: NDArray[np.floating]
) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
    """dmr_cio_91206-91236.fits contains data from day 206 of 1991 to day 236 of 1991)."""

    N_EARTH_INTERP_STEPS = 100

    dmr_day = dirbe_utils.dirbe_day_to_dmr_day(day)

    for file in dirbe_utils.DIRBE_POS_FILES:
        start_date, stop_date = dirbe_utils.get_dmrfile_datetimes(file)
        if start_date <= dmr_day <= stop_date:
            dmr_file = file
            break
    else:
        raise FileNotFoundError(
            f"could not find file containing sat pos for {dmr_day=}"
        )

    dmr_times = dirbe_utils.get_dmrfile_mjd_times(dmr_file)
    dmr_positions = dirbe_utils.get_dmrfile_positions(dmr_file)

    interpolator_sat_x = interpolate.interp1d(
        dmr_times, dmr_positions[:, 0], fill_value="extrapolate"
    )
    interpolator_sat_y = interpolate.interp1d(
        dmr_times, dmr_positions[:, 1], fill_value="extrapolate"
    )
    interpolator_sat_z = interpolate.interp1d(
        dmr_times, dmr_positions[:, 2], fill_value="extrapolate"
    )

    pos_x = interpolator_sat_x(dirbe_times)
    pos_y = interpolator_sat_y(dirbe_times)
    pos_z = interpolator_sat_z(dirbe_times)

    celestial_sat_pos = u.Quantity([pos_x, pos_y, pos_z], u.m).to(u.AU).value

    rotator = hp.Rotator(coord=["C", "E"])
    geocentric_ecl_sat_pos = u.Quantity(rotator(celestial_sat_pos), u.AU).transpose()

    interp_time_range = np.linspace(
        dirbe_times[0], dirbe_times[-1], N_EARTH_INTERP_STEPS
    )

    earth_pos = get_body("earth", Time(interp_time_range, format="mjd")).transform_to(
        HeliocentricMeanEcliptic
    )
    earth_pos = earth_pos.cartesian.xyz.to(u.AU).transpose()

    interpolator_earth_x = interpolate.interp1d(interp_time_range, earth_pos[:, 0])
    interpolator_earth_y = interpolate.interp1d(interp_time_range, earth_pos[:, 1])
    interpolator_earth_z = interpolate.interp1d(interp_time_range, earth_pos[:, 2])

    earth_pos_x = interpolator_earth_x(dirbe_times)
    earth_pos_y = interpolator_earth_y(dirbe_times)
    earth_pos_z = interpolator_earth_z(dirbe_times)

    earth_helio_pos = u.Quantity(
        [earth_pos_x, earth_pos_y, earth_pos_z], u.AU
    ).transpose()

    ecl_sat_pos = earth_helio_pos + geocentric_ecl_sat_pos

    return ecl_sat_pos.value, earth_helio_pos.value


def get_bandpass(band: int) -> tuple[u.Quantity[u.micron], NDArray[np.float64]]:
    bandpass_file = BANDPASS_PATH / f"DIRBE_{band:02}_bandpass.dat"
    bandpass = np.loadtxt(bandpass_file, unpack=True)

    non_zero = np.nonzero(bandpass[1])
    bandpass = bandpass[:, non_zero[0]]
    freqs, weights = bandpass
    freqs *= u.micron

    return freqs, weights


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


def main() -> None:
    version = 11
    smooth_pixels = True
    nside_out = 256
    print(
        f"Writing tods: {smooth_pixels=}, {nside_out=}, {version=}, {TEMP_OUTPUT_PATH=}"
    )
    write_dirbe_commmander_tods(
        output_path=TEMP_OUTPUT_PATH,
        version=version,
        smooth_pixels=smooth_pixels,
        nside_out=nside_out,
    )


if __name__ == "__main__":
    main()
