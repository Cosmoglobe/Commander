"""
dir=chains_test_$d
Final version of the commander h5 file script.

This code generates 10 h5 files, one for each dirbe band, each containing 285 pid (chunks), one
for each cio file which contains observations for 1 yday. 

Notes:
- Time gaps are filled with BAD_DATA sentinels to keep the tods continuous.
- No polarization data is included. Psi is currently as of 02-05-2023 not included in the cio files,
  but this is on the todo list.

Run as:
    OMP_NUM_THREADS=1 python write_tods_final.py 

MAKE SURE TO CHANGE VERSION NUMBER AND OUTPUT PATH BEFORE RUNNING TO NOT ACCIDENTALLY OVERWRITE EXISTING FILES (OR SET overwrite=False).
"""

from __future__ import annotations

from dataclasses import dataclass
import itertools
from pathlib import Path
from datetime import timedelta
import multiprocessing

import time
from astropy.io import fits
import astropy.units as u
import healpy as hp
from matplotlib import pyplot as plt
import quadcube
import numpy as np
import dirbe_utils
from scipy.interpolate import interp1d
from astropy.time import Time, TimeDelta
from cosmoglobe.tod_tools import TODLoader

# Path objects
DIRBE_DATA_PATH = Path("/mn/stornext/d5//data/duncanwa/DIRBE/hdf_files/")
BANDPASS_PATH = Path("/mn/stornext/d5/data/metins/dirbe/data/")
CIO_PATH = Path("/mn/stornext/d16/cmbco/ola/dirbe/cio")

# system constants
N_PROC = multiprocessing.cpu_count() // 4

ROTATOR = hp.Rotator(coord=["E", "G"])
YDAYS = np.concatenate([np.arange(89345, 89366), np.arange(90001, 90265)])
START_TIME = Time("1981-01-01", format="isot", scale="utc")

# CIO constants
N_CIO_FILES = 285
BAD_DATA_SENTINEL = -16375
TSCAL = 2e-15
SAMP_RATE = 1 / 8
SAMP_RATE_DAYS = SAMP_RATE / (24 * 3600)

BEAM_DATA = dirbe_utils.get_beam_data()

FLAG_BITS: dict[str, int] = {
    # Radiation Zone flags
    "north_van_allen_belt": 0,
    "south_van_allen_belt": 1,
    "south_atlantic_anomaly": 2,
    # Excess noise flag
    "excess_noise": 3,
    # Bad data flag
    "bad_data": 4,
    # Planet flags
    "moon": 5,
    "mercury": 6,
    "venus": 7,
    "mars": 8,
    "jupiter": 9,
    "saturn": 10,
    "uranus": 11,
    "neptune": 12,
    "comets": 13,
    "asteroids": 14,
    # Orbit and attitude flag (each observation has either of each pair turned on)
    "non_definitive_attitude": 15,
    "coarse_attitude": 16,
}

'''
The default DIRBE flags used for Commander analysis are in the defaults directory. For example,
for band 1 the default value is BAND_TOD_FLAG&&& = 2047, which includes all the radiation zone flags,
the excess noise flag, the bad data flag, and the planet flags. Currently it seems that the attitude flags
are not necessary, but they are included now for completeness.

Only a few bands are sensitive to the infrared radiation from comets and asteroids, specifically bands 4 and 5. For these,
the flags should include the comet and asteroid flags. When summed together, these flags should be 32767 = 2**15 - 1.

Comet flags and asteroid flags currently are set by hand, since their FWHM are not documented anywhere specifically, although
Arendt 2014 does specify that the comet tails are up to 15 degrees long when they are within 2 AU of the sun.
'''


@dataclass
class YdayData:
    tods: dict[str, np.ndarray]
    pixels: dict[str, np.ndarray]
    psis: dict[str, np.ndarray]
    flags: dict[str, np.ndarray]
    time_start: float
    time_stop: float
    sat_pos_start: np.ndarray
    sat_pos_stop: np.ndarray
    earth_pos_start: np.ndarray
    earth_pos_stop: np.ndarray


@dataclass
class CIO:
    tods: list[np.ndarray]
    pixels: list[np.ndarray]
    psis: list[np.ndarray]
    flags: list[np.ndarray]
    time_start: list[float]
    time_stop: list[float]
    sat_pos_start: list[np.ndarray]
    sat_pos_stop: list[np.ndarray]
    earth_pos_start: list[np.ndarray]
    earth_pos_stop: list[np.ndarray]


def get_cios(yday_data: list[YdayData]) -> dict[str, CIO]:
    return {
        f"{band:02}": CIO(
            tods=[yday.tods[f"{band:02}"] for yday in yday_data],
            pixels=[yday.pixels[f"{band:02}"] for yday in yday_data],
            psis=[yday.psis[f"{band:02}"] for yday in yday_data],
            flags=[yday.flags[f"{band:02}"] for yday in yday_data],
            time_start=[yday.time_start for yday in yday_data],
            time_stop=[yday.time_stop for yday in yday_data],
            sat_pos_start=[yday.sat_pos_start for yday in yday_data],
            sat_pos_stop=[yday.sat_pos_stop for yday in yday_data],
            earth_pos_start=[yday.earth_pos_start for yday in yday_data],
            earth_pos_stop=[yday.earth_pos_stop for yday in yday_data],
        )
        for band in dirbe_utils.BANDS
    }


def get_yday_data(
    files: range, nside_out: int, planet_time_delta: timedelta, color_corr: bool
) -> list[YdayData]:
    comet_interps, asteroid_interps = dirbe_utils.get_smallbody_interps(planet_time_delta)
    planet_interps = dirbe_utils.get_planet_interps(planet_time_delta)

    with multiprocessing.Pool(processes=N_PROC) as pool:
        proc_chunks = [
            pool.apply_async(
                get_yday_cio_data,
                args=(file_number, nside_out, planet_interps, comet_interps,
                    asteroid_interps, color_corr),
            )
            for file_number in files
        ]
        return [result.get() for result in proc_chunks if result]


def get_yday_cio_data(
    file_number: int,
    nside_out: int,
    planet_interps: dict[str, dict[str, interp1d]],
    comet_interps: dict[str, dict[str, interp1d]],
    asteroid_interps: dict[str, dict[str, interp1d]],
    color_corr: bool,
) -> YdayData:
    """Function which extracts and reorders the CIO data from one day CIO file."""

    yday = YDAYS[file_number]
    file = CIO_PATH / f"DIRBE_CIO_P3B_{yday:5}.FITS"

    try:
        with fits.open(file) as hdul:
            data = hdul[1].data

    except FileNotFoundError:
        print(f"File {file} not found")
        exit()

    # Extract time array and sort according to MJD time
    time = data["Time"].astype(np.float64)
    time_sorted_inds = np.argsort(time)
    time = time[time_sorted_inds]

    # Convert time to MJD
    time = (START_TIME + TimeDelta(time, format="sec", scale="tai")).mjd
    sat_pos_start, earth_pos_start= dirbe_utils.get_sat_and_earth_pos(yday, time[0])
    sat_pos_stop, earth_pos_stop = dirbe_utils.get_sat_and_earth_pos(yday, time[-1])

    # Gap filling
    # Get indexes where time difference is larger than 2 sampling rates and split time array
    # into chunks at those indexes.
    time_diff = np.diff(time)
    split_inds = np.nonzero(np.abs(time_diff) >= 2 * SAMP_RATE_DAYS)[0] + 1
    time_chunks = np.split(time, split_inds)

    # Create padding arrays which are to be appended at the end of the splitted chunks
    # Time is padded such that the full time array is continuous, bad data is padded with
    # BAD_DATA_SENTINEL, pixels are padded with 0 and flags are padded with 2**4 (bad_data bit).
    base_padding = [
        np.arange(chunk1[-1], chunk2[0], SAMP_RATE_DAYS)
        for chunk1, chunk2 in itertools.pairwise(time_chunks)
    ]

    bad_data_padding = padd_vals(base_padding, BAD_DATA_SENTINEL)
    pix_padding = padd_vals(base_padding, 0, dtype=np.int64)
    flag_padding = padd_vals(base_padding, 2 ** FLAG_BITS["bad_data"], dtype=np.int32)
    psi_padding = padd_vals(base_padding, 0, dtype=np.int64)

    # Get cio flags
    rad_zone_flags = data["RadZone"].astype(np.int32)[time_sorted_inds]
    xs_noise_flags = data["XSNoise"].astype(np.int32)[time_sorted_inds]
    oa_flags = data["OA_Flags"].astype(np.int32)[time_sorted_inds]

    common_flags = np.zeros_like(rad_zone_flags, dtype=np.int32)
    non_zero_rad_zone_inds = rad_zone_flags > 0
    common_flags[non_zero_rad_zone_inds > 0] += (
        2 ** rad_zone_flags[non_zero_rad_zone_inds]
    )

    common_flags += get_oa_flags(oa_flags, yday)
    # Get nside highres pixels
    pix9 = data["Pixel_no"].astype(np.int64)
    pix9_sub = data["PSubPos"].astype(np.int64)
    pix9_sub_sub = data["PSbSbPos"].astype(np.int64)
    pix15 = 4096 * pix9 + 16 * pix9_sub + pix9_sub_sub

    attack_vecs = data["AttackV"].astype(np.float64)
    attack_vecs[attack_vecs * TSCAL >= TSCAL] += 0.5 * TSCAL
    attack_vecs[attack_vecs * TSCAL <= -TSCAL] -= 0.5 * TSCAL

    natv = attack_vecs / np.expand_dims(np.linalg.norm(attack_vecs, axis=1), axis=1)

    # Get unit vectors from pixel 15 centers (uses quadcube package)
    los = quadcube.pix2vec(pix15)

    # Extract tods, and modify pointing vectors per detector according to beam data
    tods: dict[str, np.ndarray] = {}
    pixels: dict[str, np.ndarray] = {}
    psis: dict[str, np.ndarray] = {}
    flags: dict[str, np.ndarray] = {}
    for band in dirbe_utils.BANDS:
        band_label = f"{band:02}"
        beam_det_label = f"{band:01}A" if band <= 3 else f"{band:01}"
        cio_band_label = f"{band:01}A" if band <= 3 else f"{band:02}"

        # Get band specific flags
        band_bit = dirbe_utils.band_to_bit(band)

        # Pick correct beam file from yday (kinda hacky)
        for r in BEAM_DATA[beam_det_label].beam:
            if yday in r:
                isco, xsco = BEAM_DATA[beam_det_label].beam[r]
                break
        else:
            raise ValueError(f"No beam data for {yday} and {cio_band_label}")

        # Apply attack vector corrections to pointing
        new_los = np.zeros_like(los)
        new_los[0] = (
            los[0]
            + isco * natv[:, 0]
            + xsco * (los[1] * natv[:, 2] - los[2] * natv[:, 1])
        )
        new_los[1] = (
            los[1]
            + isco * natv[:, 1]
            + xsco * (los[2] * natv[:, 0] - los[0] * natv[:, 2])
        )
        new_los[2] = (
            los[2]
            + isco * natv[:, 2]
            + xsco * (los[0] * natv[:, 1] - los[1] * natv[:, 0])
        )

        # Rotate ecliptic vectors to galactic (commander convention) and apply time ordered sorting
        unit_vectors = new_los[:, time_sorted_inds]
        unit_vectors_gal = ROTATOR(unit_vectors)

        # Convert unit vectors to requested nside resolution healpix pixels
        pix = hp.vec2pix(nside_out, *unit_vectors_gal)
        lon, lat = hp.pix2ang(nside_out, pix, lonlat=True)
        v_x, v_y, v_z = unit_vectors_gal

        tod = data[f"Phot{cio_band_label}"].astype(np.float64)[time_sorted_inds]

        # Compute polarization angle
        lonrad = np.deg2rad(lon)
        latrad = np.deg2rad(lat)
        colat = np.pi/2 - latrad
        att_x = natv[:,0]*np.cos(lonrad) + natv[:,1]*np.sin(lonrad)
        att_y = natv[:,1]*np.cos(lonrad) - natv[:,0]*np.sin(lonrad)

        xp = att_x*np.cos(colat) - natv[:,2]*np.sin(colat)
        psi = np.arctan2(att_y, xp)
        psi[psi < 0] += 2*np.pi
        psi = psi % (2*np.pi)


        # Remove the iras convention color correction
        iras_color_corr_factor = dirbe_utils.get_iras_factor(band)
        flag = common_flags.copy()
        # Get planet flags
        for body, radius in dirbe_utils.PLANET_RADII.items():
            planet_lon = planet_interps[body]["lon"](time)
            planet_lat = planet_interps[body]["lat"](time)

            ang_dist = hp.rotator.angdist(
                np.array([lon, lat]),
                np.array([planet_lon, planet_lat]),
                lonlat=True,
            )

            if (body == 'moon') & ((band == 4) | (band == 5)):
                r = 15
            else:
                r = radius
            planet_indices = ang_dist <= np.deg2rad(r)
            flag[planet_indices] += 2 ** FLAG_BITS[body]


        # Get comet flags
        for body, radius in dirbe_utils.COMET_RADII.items():
            comet_x = comet_interps[body]["x"](time)
            comet_y = comet_interps[body]["y"](time)
            comet_z = comet_interps[body]["z"](time)
            comet_dist = comet_interps[body]["dist"](time)
            if comet_dist.min() < 2:
                r = radius
            else:
                r = 1

            ang_dist = hp.rotator.angdist(
                np.array([v_x, v_y, v_z]),
                np.array([comet_x, comet_y, comet_z])
            )
            

            comet_indices = ang_dist <= np.deg2rad(r)
            flag[comet_indices] += 2 ** FLAG_BITS["comets"]

        for body, radius in dirbe_utils.ASTEROID_RADII.items():
            aster_x = asteroid_interps[body]["x"](time)
            aster_y = asteroid_interps[body]["y"](time)
            aster_z = asteroid_interps[body]["z"](time)

            ang_dist = hp.rotator.angdist(
                np.array([v_x, v_y, v_z]),
                np.array([aster_x, aster_y, aster_z])
            )

            aster_indices = ang_dist <= np.deg2rad(radius)
            flag[aster_indices] += 2 ** FLAG_BITS["asteroids"]

        # Find which flags correspond to the detector and get the inds of those flags
        xs_noise_inds = (xs_noise_flags & band_bit) > 0
        flag[xs_noise_inds] += 2 ** FLAG_BITS["excess_noise"]

        # Get bad data flags
        bad_data_inds = tod <= BAD_DATA_SENTINEL
        flag[bad_data_inds] += 2 ** FLAG_BITS["bad_data"]

        # padd tods, pix and flags to remove gaps in data
        pixels[band_label] = padd_array_gaps(
            np.split(pix, split_inds), padding=pix_padding
        )

        tods[band_label] = padd_array_gaps(
            np.split(tod * iras_color_corr_factor if color_corr else tod, split_inds),
            padding=bad_data_padding,
        )

        psis[band_label] = padd_array_gaps(
            np.split(psi, split_inds), padding=psi_padding
        )

        flags[band_label] = padd_array_gaps(
            np.split(flag, split_inds), padding=flag_padding
        )

    return YdayData(
        tods, 
        pixels, 
        psis,
        flags, 
        time_start=time[0], 
        time_stop=time[-1], 
        sat_pos_start=sat_pos_start, 
        sat_pos_stop=sat_pos_stop, 
        earth_pos_start=earth_pos_start,
        earth_pos_stop=earth_pos_stop,
    )


def padd_array_gaps(splits: list[np.ndarray], padding: list[np.ndarray]) -> np.ndarray:
    return np.concatenate([np.append(s, p) for s, p in zip(splits, padding)])


def padd_vals(
    base_padding: list[np.ndarray], val: int | float, dtype: np.dtype = float
) -> list[np.ndarray]:
    return [np.full_like(padding, val, dtype=dtype) for padding in base_padding]


def get_oa_flags(oa_flags: np.ndarray, yday: int) -> np.ndarray:
    """
    Function which converts the OA flags to the new flag system. The new flag system is as follows:
    - The OA flags are split into 7 bits, each bit corresponding to a different flag.
    From the DIRBE Explanatory supplement, these are
    - Bit 0: 1 if non-definitive attitude, 0 if definitive attitude
    - Bit 1: 1 if coarse attitude, 0 if fine attitude
    - Bit 2: 1 if merged attitude, 0 if external UAX attitude
    - Bit 3: 1 if spacecraft slewing, 0 if spacecraft not slewing
    - Bit 4: 1 if special pointing, 0 if normal pointing
    - Bit 5: 1 when spacecraft in ascending portion of COBE orbit, 
                i.e., moving up toward North Ecliptic Pole,
             0 when spacecraft in descending portion of COBE orbit,
                i.e., moving down toward South Ecliptic Pole
    - Bit 6: 1 if leading LOS, 0 if trailing LOS
    - Bit 7: Not used

    Currently, only the non-definitive attitude is considered a gamechanger, so flags with value of 15 (2**15 = 32768)
    need to be removed. All others are set to 2**16 = 65536

    Default assumptions si taht the attitude solutions are the highest values. 
    """
    new_bits = iter(
        [
            bit
            for bit in FLAG_BITS.values()
            if bit >= FLAG_BITS["non_definitive_attitude"]
        ]
    )
    flags = np.zeros_like(oa_flags)
    for cio_bit, bit1 in zip(range(7), new_bits):
        inds = (oa_flags & 2**cio_bit) > 0
        flags[inds] += 2**bit1
    return flags


def get_flag_sum(flags: list[str]) -> int:

    if flags:
        return sum(2 ** FLAG_BITS[flag] for flag in flags)
    return 0


def write_band(
    comm_tod: TODLoader, cio: CIO, filename: str, band: int, nside_out: int, n_pids: int
) -> None:
    COMMON_GROUP = "/common"
    HUFFMAN_COMPRESSION = ["huffman", {"dictNum": 1}]

    det_str = f"{band:02}_A"
    comm_tod.init_file(freq=filename, od="", mode="w")

    comm_tod.add_field(COMMON_GROUP + "/fsamp", 1 / SAMP_RATE)
    comm_tod.add_field(COMMON_GROUP + "/nside", [nside_out])
    comm_tod.add_field(COMMON_GROUP + "/det", np.string_(det_str + ","))

    comm_tod.add_field(COMMON_GROUP + "/polang", [0])
    comm_tod.add_attribute(COMMON_GROUP + "/polang", "index", det_str)

    comm_tod.add_field(COMMON_GROUP + "/mbang", [0])
    comm_tod.add_attribute(COMMON_GROUP + "/mbang", "index", det_str)

    for pid in range(n_pids):
        pid_label = f"{pid+1:06}"
        pid_common_group = pid_label + "/common"
        pid_det_group = f"{pid_label}/{det_str}"

        comm_tod.add_field(pid_common_group + "/time", [cio.time_start[pid], 0, 0])
        comm_tod.add_attribute(pid_common_group + "/time", "index", "MJD, OBT, SCET")

        comm_tod.add_field(pid_common_group + "/time_end", [cio.time_stop[pid], 0, 0])
        comm_tod.add_attribute(pid_common_group + "/time_end", "index", "MJD, OBT, SCET")

        comm_tod.add_field(pid_common_group + "/ntod", [len(cio.tods[pid])])

        comm_tod.add_field(pid_common_group + "/satpos", cio.sat_pos_start[pid])
        comm_tod.add_field(pid_common_group + "/satpos_end", cio.sat_pos_stop[pid])

        comm_tod.add_field(pid_common_group + "/earthpos", cio.earth_pos_start[pid])
        comm_tod.add_field(pid_common_group + "/earthpos_end", cio.earth_pos_stop[pid])

        comm_tod.add_attribute(pid_common_group + "/satpos", "index", "X, Y, Z")
        comm_tod.add_attribute(
            pid_common_group + "/satpos", "coords", "heliocentric-ecliptic"
        )

        comm_tod.add_field(pid_det_group + "/flag", cio.flags[pid], HUFFMAN_COMPRESSION)

        comm_tod.add_field(pid_det_group + "/tod", cio.tods[pid])
        comm_tod.add_field(pid_det_group + "/pix", cio.pixels[pid], HUFFMAN_COMPRESSION)

        # TODO: Get correct polarization angle (detector angle)
        psi_digitize_compression = [
            "digitize",
            {"min": 0, "max": 2 * np.pi, "nbins": 64},
        ]
        comm_tod.add_field(
            pid_det_group + "/psi",
            cio.psis[pid],
            [psi_digitize_compression, HUFFMAN_COMPRESSION],
        )

        comm_tod.add_field(pid_det_group + "/outP", np.zeros((2, 1)))

        const_scalars = dirbe_utils.get_const_scalars(band)
        comm_tod.add_field(pid_det_group + "/scalars", const_scalars)
        comm_tod.add_attribute(
            pid_det_group + "/scalars", "index", "gain, sigma0, fknee, alpha"
        )

        comm_tod.finalize_chunk(pid + 1)

    comm_tod.finalize_file()


def write_to_commander_tods(
    cios: dict[str, CIO],
    nside_out: int,
    version: int,
    out_path: Path,
    overwrite: bool = False,
) -> None:
    manager = multiprocessing.Manager()
    pool = multiprocessing.Pool(processes=N_PROC)

    multiprocessor_manager_dicts = {}
    for band in dirbe_utils.BANDS:
        # name = f"DIRBE_{band:02}_nside{nside_out:03}_zodi_only"
        name = f"DIRBE_{band:02}_nside{nside_out:03}_V{version:02}"
        multiprocessor_manager_dicts[name] = manager.dict()

    filenames = list(multiprocessor_manager_dicts.keys())
    comm_tod = TODLoader(
        out_path, "", version, multiprocessor_manager_dicts, overwrite=overwrite
    )

    n_pids = 0
    for cio in cios.values():
        n_pids = len(cio.time_start)
        break
    if n_pids == 0:
        raise ValueError("No CIOs found")

    # Currently, writing commander h5 files uses a pretty unoptimal interface which is difficult to
    # parallelize without corrupting the written files. This is the go to way to concurrently write files
    # which is pretty ugly.
    x = [
        [
            pool.apply_async(
                write_band,
                args=(
                    comm_tod,
                    cios[f"{band:02}"],
                    filenames[band - 1],
                    band,
                    nside_out,
                    n_pids,
                ),
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


def main() -> None:
    time_delta = timedelta(hours=1)
    files = range(N_CIO_FILES)
    nside_out = 512

    start_time = time.perf_counter()
    color_corr = False
    version = 21

    print(f"{'Writing DIRBE h5 files':=^50}")
    print(f"{version=}, {nside_out=}")
    print(f"reading and processing cios for {len(files)} ydays...")
    yday_data = get_yday_data(
        files, nside_out=nside_out, planet_time_delta=time_delta, color_corr=color_corr
    )
    cios = get_cios(yday_data)
    cio_time = time.perf_counter() - start_time
    print("done")
    print(
        f"time spent reading in and preprocessing cios: {(cio_time/60):2.2f} minutes\n"
    )

    print("writing cios to h5 files...")
    write_to_commander_tods(
        cios,
        nside_out=nside_out,
        version=version,
        out_path=DIRBE_DATA_PATH,
        overwrite=True,
    )
    h5_time = time.perf_counter() - start_time
    print("done")
    print(f"time spent writing to h5: {(h5_time/60):2.2f} minutes\n")
    print(f"total time: {((h5_time + cio_time)/60):2.2f} minutes")
    print(f"{'':=^50}")
    exit()
    # print(cio.time.shape)
    # print(cio.tod["04"].shape)
    # import matplotlib.pyplot as plt
    # plt.plot(cio.time, cio.tod["04"])
    # plt.show()

    print("Binning map..")
    binned_map = np.zeros(hp.nside2npix(nside_out))
    total_hits = np.zeros_like(binned_map)
    cio = cios["07"]
    pix = np.concatenate([pix for pix in cio.pixels])
    tod = np.concatenate([tod for tod in cio.tods])
    flags = np.concatenate([flag for flag in cio.flags])

    flag_bit_sum = get_flag_sum(
        [
            "north_van_allen_belt",
            "south_van_allen_belt",
            "south_atlantic_anomaly",
            "excess_noise",
            "bad_data",
            "moon",
            "mercury",
            "venus",
            "mars",
            "jupiter",
            "saturn",
            "uranus",
            "neptune",
            "non_definitive_attitude",
            "coarse_attitude",
            "special_pointing",
        ]
    )
    print(f"flag bit sum: {flag_bit_sum}")

    if flag_bit_sum > 0:
        condition = np.bitwise_and(flags, flag_bit_sum) == 0
        pix = pix[condition]
        tod = tod[condition]

    bin_count = np.bincount(pix, weights=tod, minlength=len(binned_map))

    binned_map[: len(bin_count)] += bin_count
    unique_pix, counts = np.unique(pix, return_counts=True)
    total_hits[unique_pix] += counts
    non_zero_inds = total_hits > 0
    binned_map[non_zero_inds] /= total_hits[non_zero_inds]
    binned_map[binned_map <= 0] = 0
    # hp.write_map("flagged.fits", binned_map, overwrite=True)
    hp.mollview(binned_map, norm="hist")
    plt.show()
    # binned_map *= iras_factor
    # hp.mollview(binned_map, min=0, max=5)
    # plt.show()


if __name__ == "__main__":
    main()
