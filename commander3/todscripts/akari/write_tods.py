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
import akari_utils
from scipy.interpolate import interp1d
from astropy.time import Time, TimeDelta
from cosmoglobe.tod_tools import TODLoader
import pickle
import zodipy

# zodi_model = zodipy.Zodipy(extrapolate=True)
# Path objects
AKARI_DATA_PATH = Path("/mn/stornext/d5/data/duncanwa/akari/data")
HDF5_PATH = Path("/mn/stornext/d16/cmbco/bp/gustavbe/master/dirbe_hdf5_files/")
CIO_PATH = Path("/mn/stornext/u3/hke/data_d5/akari/TOD")

# system constants
N_PROC = multiprocessing.cpu_count()

ROTATOR = hp.Rotator(coord=["E", "G"])
YDAYS = np.concatenate([np.arange(89345, 89366), np.arange(90001, 90265)])
START_TIME = Time("1981-01-01", format="isot", scale="utc")

# CIO constants
N_CIO_FILES = 285
BAD_DATA_SENTINEL = -16375
TSCAL = 2e-15
SAMP_RATE = 1 / 8
SAMP_RATE_DAYS = SAMP_RATE / (24 * 3600)

#BEAM_DATA = akari_utils.get_beam_data()


# flags exist for each detector.

# bad_frame, undef_anom_Frame, blank, in_saa, near_moon, untrusted_frame,
# reserve. These are for the entire frame, all elements


# See page 63 of https://www.ir.isas.jaxa.jp/AKARI/Observation/support/FIS/IDUM/FIS_IDUM_1.3.pdf

# flag is for each frame
# bit_flag is for each frame, 32 flags per detector
# status 

# calala - means "calibration lamp"

# The flags set in general are
# dead      gpl_type1 mtgl_type1
# saturate  gpl_type2 mtgl_type2
# reset     gpl_type3 mtgl_type2
# rstanom   gpl_type4 mtgl_type3
# no_diff   gpl_tail


# dead pixels don't need to be included

# Integer ADU values are included.

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
    # Orbit and attitude flag (each observation has either of each pair turned on)
    "non_definitive_attitude": 13,
    "definite_attitude": 14,
    "course_attitude": 15,
    "fine_attitude": 16,
    "merged_attitude": 17,
    "external_uax_attitude": 18,
    "space_craft_slewing": 19,
    "space_craft_not_slewing": 20,
    "special_pointing": 21,
    "normal_pointing": 22,
    "space_craft_ascending": 23,
    "space_craft_descending": 24,
    "leading_los": 25,
    "trailing_los": 26,
}


@dataclass
class YdayData:
    tods: dict[str, np.ndarray]
    pixels: dict[str, np.ndarray]
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
    flags: list[np.ndarray]
    time_start: list[float]
    time_stop: list[float]
    sat_pos_start: list[np.ndarray]
    sat_pos_stop: list[np.ndarray]
    earth_pos_start: list[np.ndarray]
    earth_pos_stop: list[np.ndarray]


def get_cios(yday_data: list[YdayData]) -> dict[str, CIO]:
    output_dict = {}
    for band in akari_utils.DETECTORS:
        for yday in yday_data:
            if f'Akari_{band}' in yday.tods.keys():
                output_dict[f'Akari_{band}'] =  CIO(tods=[yday.tods[f"Akari_{band}"]],
                     pixels=[yday.pixels[f"Akari_{band}"]],
                     flags=[yday.flags[f"Akari_{band}"]],
                     time_start=[yday.time_start],
                     time_stop=[yday.time_stop],
                     sat_pos_start=[yday.sat_pos_start],
                     sat_pos_stop=[yday.sat_pos_stop],
                     earth_pos_start=[yday.earth_pos_start],
                     earth_pos_stop=[yday.earth_pos_stop])
    return output_dict


def get_yday_data(
        tod: NDArray, lon: NDArray, lat: NDArray, nside_out: int, planet_time_delta: timedelta, color_corr: bool
) -> list[YdayData]:

    bands = ['WIDE-S']*60 + ['N60']*40
    inds = np.concatenate((np.arange(60), np.arange(40))) + 1
    data_ind = np.arange(100)

    with multiprocessing.Pool(processes=N_PROC) as pool:
        proc_chunks = [
            pool.apply_async(
                get_yday_cio_data,
                args=(tod[di], lon[di], lat[di], band, ind, nside_out, color_corr),
            )
            for band, ind, di  in zip(bands, inds, data_ind)
        ]
        return [result.get() for result in proc_chunks if result]


def get_yday_cio_data(
    tod: NDArray,
    lon: NDArray,
    lat: NDArray,
    band: str,
    ind: int,
    nside_out: int,
    #planet_interps: dict[str, dict[str, interp1d]],
    color_corr: bool,
) -> YdayData:
    """Function which extracts and reorders the CIO data from one day CIO file."""


    #yday = YDAYS[file_number]

    time = 53826

    # Convert time to MJD
    time = (START_TIME + TimeDelta(time, format="sec", scale="tai")).mjd
    sat_pos_start, earth_pos_start = 0,0
    sat_pos_stop, earth_pos_stop = 0,0


    # Extract tods, and modify pointing vectors per detector according to beam data
    tods: dict[str, np.ndarray] = {}
    pixels: dict[str, np.ndarray] = {}
    flags: dict[str, np.ndarray] = {}
    band_label = f'Akari_{band}_{ind:02}'


    # padd tods, pix and flags to remove gaps in data
    pixels[band_label] = hp.ang2pix(nside_out, lon, lat, lonlat=True)

    tods[band_label] = tod

    flags[band_label] = pixels[band_label]*0

    return YdayData(
        tods, 
        pixels, 
        flags, 
        time_start=time, 
        time_stop=time, 
        sat_pos_start=sat_pos_start, 
        sat_pos_stop=sat_pos_stop, 
        earth_pos_start=earth_pos_start,
        earth_pos_stop=earth_pos_stop,
    )


def padd_array_gaps(splits: list[np.ndarray], padding: list[np.ndarray]) -> np.ndarray:
    return np.concatenate([np.append(s, p) for s, p in zip(splits, padding)])





def write_band(
        comm_tod: TODLoader, cio: CIO, filename: str, band: str, ndet: int, nside_out: int, n_pids: int
) -> None:
    COMMON_GROUP = "/common"
    HUFFMAN_COMPRESSION = ["huffman", {"dictNum": 1}]

    det_str = ""
    for i in range(ndet):
        det_str = det_str + f'{band}_{i+1:02}, '
    comm_tod.init_file(freq=filename, od="", mode="w")



    for v in cio.keys():
        det1 = v
        break

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

        comm_tod.add_field(pid_common_group + "/time", [cio[det1].time_start[pid], 0, 0])
        comm_tod.add_attribute(pid_common_group + "/time", "index", "MJD, OBT, SCET")

        comm_tod.add_field(pid_common_group + "/time_end", [cio[det1].time_stop[pid], 0, 0])
        comm_tod.add_attribute(pid_common_group + "/time_end", "index", "MJD, OBT, SCET")

        comm_tod.add_field(pid_common_group + "/ntod", [len(cio[det1].tods[pid])])

        comm_tod.add_field(pid_common_group + "/satpos", cio[det1].sat_pos_start[pid])
        comm_tod.add_field(pid_common_group + "/satpos_end", cio[det1].sat_pos_stop[pid])

        comm_tod.add_field(pid_common_group + "/earthpos", cio[det1].earth_pos_start[pid])
        comm_tod.add_field(pid_common_group + "/earthpos_end", cio[det1].earth_pos_stop[pid])

        comm_tod.add_attribute(pid_common_group + "/satpos", "index", "X, Y, Z")
        comm_tod.add_attribute(
            pid_common_group + "/satpos", "coords", "heliocentric-ecliptic"
        )

        for i in range(ndet):
            det_lab = f'Akari_{band}_{i+1:02}'
            pid_det_group = f"{pid_label}/{i+1:02}"
            comm_tod.add_field(pid_det_group + "/flag", cio[det_lab].flags[pid], HUFFMAN_COMPRESSION)

            comm_tod.add_field(pid_det_group + "/tod", cio[det_lab].tods[pid])
            comm_tod.add_field(pid_det_group + "/pix", cio[det_lab].pixels[pid], HUFFMAN_COMPRESSION)

            # TODO: Get correct polarization angle (detector angle)
            psi_digitize_compression = [
                "digitize",
                {"min": 0, "max": 2 * np.pi, "nbins": 64},
            ]
            comm_tod.add_field(
                pid_det_group + "/psi",
                np.zeros_like(cio[det_lab].tods[pid]),
                [psi_digitize_compression, HUFFMAN_COMPRESSION],
            )

            comm_tod.add_field(pid_det_group + "/outP", np.zeros((2, 1)))

            const_scalars = akari_utils.get_const_scalars(band)
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
    filenames = {}
    for band in akari_utils.BANDS:
        # name = f"DIRBE_{band:02}_nside{nside_out:03}_zodi_only"
        name = f"Akari_{band}_v{version:02}"
        multiprocessor_manager_dicts[name] = manager.dict()
        filenames[band] = name

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
    x = [[]]
    for ndet, band in zip(akari_utils.NDETS[:2], akari_utils.BANDS[:2]):
        print(band)
        x[0].append(
        pool.apply_async(
            write_band,
            args=(
                comm_tod,
                cios,
                filenames[band],
                band,
                ndet,
                nside_out,
                n_pids)))

    for res1 in x:
        for res in res1:
            res.get()

    pool.close()
    pool.join()

    comm_tod.make_filelists()


def main() -> None:
    time_delta = timedelta(hours=1)
    files = range(1)
    nside_out = 8196

    start_time = time.perf_counter()
    color_corr = False
    version = 0

    print(f"{'Writing Akari h5 files':=^50}")
    print(f"{version=}, {nside_out=}")
    print(f"reading and processing cios for {len(files)} ydays...")
    datfile = CIO_PATH/'FIS_SW_20060801000000_flux.pkl'
    latfile = CIO_PATH/'FIS_SW_20060801000000_gb_lat.pkl'
    lonfile = CIO_PATH/'FIS_SW_20060801000000_gb_lon.pkl'
    with open(datfile, 'rb') as f:
        tods = pickle.load(f)
    with open(latfile, 'rb') as f:
        lats = pickle.load(f)
    with open(lonfile, 'rb') as f:
        lons = pickle.load(f)

    yday_data = get_yday_data(
        tods, lons, lats, nside_out=nside_out, planet_time_delta=time_delta, color_corr=color_corr
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
        out_path=AKARI_DATA_PATH,
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
            # "definite_attitude",
            "course_attitude",
            # "fine_attitude",
            # "merged_attitude",
            # "external_uax_attitude",
            # "space_craft_slewing",
            # "space_craft_not_slewing",
            "special_pointing",
            # "normal_pointing",
            # "space_craft_ascending",
            # "space_craft_descending",
            # "leading_los",
            # "trailing_los",
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
