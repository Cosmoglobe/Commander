"""
version 2: add smoothing of pixels to remove pixelization effects
version 3: rotate tods to galactic coordinates
version 4: add custom planet flags
version 5: larger flagging radii
version 6: per pointing time for accurate planet flags
"""

from __future__ import annotations
from multiprocessing.managers import DictProxy
import os
import random

from typing import TYPE_CHECKING
import sys
from scipy import interpolate
from datetime import datetime, timedelta

if TYPE_CHECKING:
    from ...python.commander_tools.tod_tools.commander_tod import (
        commander_tod,
    )

sys.path.insert(0, "/mn/stornext/d16/cmbco/bp/metins/Commander/commander3/python")
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
from astropy.time import Time, TimeDelta

version = 7
smooth_pixels = True
nside_out = 128

TEMP_OUTPUT_PATH = "/mn/stornext/d16/cmbco/bp/metins/dirbe/data/"
PATH_TO_HDF5_FILES = "/mn/stornext/d16/cmbco/bp/gustavbe/master/dirbe_hdf5_files/"

NUM_PROCS = 10
NUM_CHUNKS_PER_BAND = 285
DEFAULT_NSIDE = 128
NPSI = 2048
FSAMP = 8
BAD_DATA_SENTINEL = -16375
N_SMOOTHING_BOUNDARY = 100
N_CONV_BOX_POINTS = 30
NSIDE = 128

FLAG_BODIES = [
    "moon",
    "mercury",
    "venus",
    "mars",
    "jupiter",
    "saturn",
    "uranus",
    "neptune",
]


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

    planet_interp_times = Time(
        [
            (datetime(1989, 9, 11) + timedelta(hours=i)).isoformat()
            for i in range(int(425 * 24))
        ],
        format="isot",
        scale="utc",
    )
    planet_interps = get_planet_lonlats(planet_interp_times)

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
                write_detector,
                args=[
                    comm_tod,
                    detector,
                    nside_in,
                    nside_out,
                    smooth_pixels,
                    filenames[detector - 1],
                    planet_interps,
                ],
            )
            for detector in dirbe_utils.DETECTORS
        ]
    ]

    for res1 in x:
        for res in res1:
            res.get()

    pool.close()
    pool.join()

    comm_tod.make_filelists()


def write_detector(
    comm_tod: commander_tod,
    detector: int,
    nside_in: int,
    nside_out: int,
    smooth: bool,
    filename: str,
    planet_interps: dict[str, dict[str, interp1d]],
) -> None:
    """Writes a single chunk of tod to file."""

    hdf5_filename = PATH_TO_HDF5_FILES + f"Phot{detector:02}_{nside_in}.hdf5"
    COMMON_GROUP = "/common"
    HUFFMAN_COMPRESSION = ["huffman", {"dictNum": 1}]

    comm_tod.init_file(freq=filename, od="", mode="w")

    comm_tod.add_field(COMMON_GROUP + "/fsamp", FSAMP)

    # nside
    comm_tod.add_field(COMMON_GROUP + "/nside", [nside_out])

    # make detector names lookup
    if detector <= 3:
        detector_names = ",".join(
            [f"{detector:02}_{band_label}" for band_label in dirbe_utils.BANDS_LABELS]
        )
    else:
        detector_names = f"{detector:02}_{dirbe_utils.BANDS_LABELS[0]},"
    comm_tod.add_field(COMMON_GROUP + "/det", np.string_(detector_names))

    common_polang = get_detector_polang(detector)
    comm_tod.add_field(COMMON_GROUP + "/polang", common_polang)
    comm_tod.add_attribute(COMMON_GROUP + "/polang", "index", detector_names)

    common_mbang = get_mbang(detector)
    comm_tod.add_field(COMMON_GROUP + "/mbang", common_mbang)
    comm_tod.add_attribute(COMMON_GROUP + "/mbang", "index", detector_names)

    rotator = hp.Rotator(coord=["E", "G"])

    for chunk in range(1, NUM_CHUNKS_PER_BAND + 1):
        day = chunk
        print(f"processing chunk: {chunk:03}...")
        chunk_label = f"{chunk:06}"

        chunk_common_group = chunk_label + "/common"

        mjd_times = get_chunk_mjd_times(chunk_label, hdf5_filename)
        comm_tod.add_field(
            chunk_common_group + "/time", [mjd_times[0], 0, 0]
        )  # TODO: get SCET and OBT times
        comm_tod.add_attribute(chunk_common_group + "/time", "index", "MJD, OBT, SCET")

        ntod = get_ntods(chunk_label, hdf5_filename)
        if smooth:
            ntod -= int(2 * N_SMOOTHING_BOUNDARY)
        comm_tod.add_field(chunk_common_group + "/ntod", [ntod])

        sat_pos, earth_pos = get_sat_and_earth_pos(day, mjd_times)
        comm_tod.add_field(chunk_common_group + "/satpos", sat_pos[0])

        # Add earth position. Required for zodiacal light calculation.
        comm_tod.add_field(chunk_common_group + "/earthpos", earth_pos[0])

        # add metadata
        comm_tod.add_attribute(chunk_common_group + "/satpos", "index", "X, Y, Z")
        comm_tod.add_attribute(
            chunk_common_group + "/satpos", "coords", "heliocentric-ecliptic"
        )

        for band in dirbe_utils.BANDS_LABELS:
            band_chunk_group = f"{chunk_label}/{detector:02}_{band}/"

            tods = get_chunk_band_tods(chunk_label, hdf5_filename, band)
            pixels = get_chunk_band_pixels(
                chunk_label,
                hdf5_filename,
                band,
                rotator,
                nside_in=nside_in,
                nside_out=nside_out,
                smooth=smooth,
            )
            flags = get_chunk_band_flags(
                chunk_label,
                hdf5_filename,
                band,
                nside_in,
                planet_interps,
                pixels,
                tods,
            )
            psi = get_chunk_band_psi(chunk_label, hdf5_filename, band)

            if smooth:  # smoothing currently requires removing boundaries
                flags = flags[N_SMOOTHING_BOUNDARY:-N_SMOOTHING_BOUNDARY]
                tods = tods[N_SMOOTHING_BOUNDARY:-N_SMOOTHING_BOUNDARY]
                pixels = pixels[N_SMOOTHING_BOUNDARY:-N_SMOOTHING_BOUNDARY]
                psi = psi[N_SMOOTHING_BOUNDARY:-N_SMOOTHING_BOUNDARY]
            if flags.size > 0:
                comm_tod.add_field(
                    band_chunk_group + "/flag", flags, HUFFMAN_COMPRESSION
                )

            if tods.size > 0:
                comm_tod.add_field(band_chunk_group + "/tod", tods)

            if pixels.size > 0:
                comm_tod.add_field(
                    band_chunk_group + "/pix", pixels, HUFFMAN_COMPRESSION
                )

            psi_digitize_compression = [
                "digitize",
                {"min": 0, "max": 2 * np.pi, "nbins": NPSI},
            ]
            if psi.size > 0:
                comm_tod.add_field(
                    band_chunk_group + "/psi",
                    psi,
                    [psi_digitize_compression, HUFFMAN_COMPRESSION],
                )

            out_ang = get_out_ang()
            comm_tod.add_field(band_chunk_group + "/outP", out_ang)

            scalars = get_scalars(tods)
            comm_tod.add_field(band_chunk_group + "/scalars", scalars)
            comm_tod.add_attribute(
                band_chunk_group + "/scalars", "index", "gain, sigma0, fknee, alpha"
            )

            if detector > 3:
                break

        comm_tod.finalize_chunk(chunk)

        # if chunk > 1:
        #     break
    comm_tod.finalize_file()


def get_detector_polang(detector: int) -> list[float]:
    """TODO: get correct angles"""

    if detector > 3:
        return [0]
    return [0, 0, 0]


def get_mbang(detector: int) -> list[float]:
    """TODO: get correct angles"""

    if detector > 3:
        return [0]
    return [0, 0, 0]


def get_chunk_mjd_times(chunk_label: str, hdf5_filename: str) -> NDArray[np.floating]:
    """Gets dirbe mjd times from gustavs files."""

    with h5py.File(hdf5_filename, "r") as file:
        return file[f"{chunk_label}/A/time"][()]


def get_planet_lonlats(times: Time) -> dict[str, dict[str, interp1d]]:
    print("precomputing planet interpolators...")
    interpolaters = {}
    with solar_system_ephemeris.set("de432s"):
        for body_name in FLAG_BODIES:
            interpolaters[body_name] = {}
            body = get_body(body_name, times).galactic

            interpolaters[body_name]["l"] = interp1d(times.mjd, body.l.value)
            interpolaters[body_name]["b"] = interp1d(times.mjd, body.b.value)

    print("done")
    return interpolaters


def get_chunk_band_flags(
    chunk_label: str,
    hdf5_filename: str,
    band: str,
    nside_in: int,
    planet_interpolaters: dict[str, dict[str, interp1d]],
    pixels: NDArray[np.integer],
    tods: NDArray[np.floating],
) -> NDArray[np.integer]:
    """Gets dirbe flags from gustavs files."""
    flagging_radii = 5

    sentinel_indices = tods <= BAD_DATA_SENTINEL
    flags[sentinel_indices] += int(2**13) # Add sentinel flag for values less than SENTINEL_VALUE

    with h5py.File(hdf5_filename, "r") as file:
        flags = file[f"{chunk_label}/{band}/flag"][()]
        time = file[f"{chunk_label}/{band}/time"][()]


        for bit_idx, body in enumerate(FLAG_BODIES, start=14):
            lon, lat = hp.pix2ang(nside_in, pixels, lonlat=True)
            planet_lon = planet_interpolaters[body]["l"](time)
            planet_lat = planet_interpolaters[body]["b"](time)
            ang_dist = hp.rotator.angdist(
                np.asarray([lon, lat]),
                np.asarray([planet_lon, planet_lat]),
                lonlat=True,
            )

            planet_indices = ang_dist <= np.deg2rad(flagging_radii)
            flags[planet_indices] += int(2 ** (bit_idx))

        return flags


def get_chunk_band_tods(
    chunk_label: str, hdf5_filename: str, band: str
) -> NDArray[np.floating]:
    """Gets dirbe tods from gustavs files."""

    with h5py.File(hdf5_filename, "r") as file:
        return file[f"{chunk_label}/{band}/tod"][()]


def get_ntods(chunk_label: str, hdf5_filename: str) -> int:
    """Gets dirbe tods from gustavs files."""

    with h5py.File(hdf5_filename, "r") as file:
        return len(file[f"{chunk_label}/A/tod"])


def get_chunk_band_psi(
    chunk_label: str, hdf5_filename: str, band: str
) -> NDArray[np.floating]:
    """Gets dirbe polang from gustavs files. TODO: rotate to galactic."""

    with h5py.File(hdf5_filename, "r") as file:
        return file[f"{chunk_label}/{band}/polang"][()]


def get_chunk_band_pixels(
    chunk_label: str,
    hdf5_filename: str,
    band: str,
    rotator: hp.Rotator,
    nside_in: int,
    nside_out: int,
    smooth: bool = False,
) -> NDArray[np.integer]:
    """Gets dirbe pixels from gustavs files."""

    with h5py.File(hdf5_filename, "r") as file:
        pixels_ecl = file[f"{chunk_label}/{band}/pix"][()]

    unit_vectors_ecl = np.asarray(hp.pix2vec(nside_in, pixels_ecl))
    if smooth:
        unit_vectors_ecl = np.array(
            [
                np.convolve(
                    unit_vector,
                    np.ones(N_CONV_BOX_POINTS) / N_CONV_BOX_POINTS,
                    mode="same",
                )
                for unit_vector in unit_vectors_ecl
            ]
        )

    unit_vectors_gal = rotator(unit_vectors_ecl)
    return hp.vec2pix(nside_out, *unit_vectors_gal)


def get_out_ang() -> NDArray[np.floating]:
    """TODO: get correct out ang maybe..?"""

    return np.zeros((2, 1))


def get_scalars(tods: NDArray[np.floating]) -> NDArray[np.floating]:
    """TODO: get correct gain, sigma0, fknee and alpha."""

    TEMP_GAIN = 1
    TEMP_ALPHA = -1
    fknee = 1 / (10 * 60)

    sigma0 = np.diff(tods).std() / (2**0.5)

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


def main() -> None:
    print(f"Writing tods: {smooth_pixels=}, {nside_out=}")
    write_dirbe_commmander_tods(
        output_path=TEMP_OUTPUT_PATH,
        version=version,
        smooth_pixels=smooth_pixels,
        nside_out=nside_out,
    )


if __name__ == "__main__":
    main()
