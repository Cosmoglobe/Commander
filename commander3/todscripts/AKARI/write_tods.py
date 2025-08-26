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
from astropy.coordinates import SkyCoord
import healpy as hp
from matplotlib import pyplot as plt
import quadcube
import numpy as np
import akari_utils
from scipy.interpolate import interp1d
from astropy.time import Time, TimeDelta
from cosmoglobe.tod_tools import TODLoader
import pickle

#from cProfile import Profile
#from pstats import SortKey, Stats

# Path objects
AKARI_DATA_PATH = Path("/mn/stornext/d5/data/duncanwa/akari/data")
AKARI_DATA_PATH = Path("/mn/stornext/d5/data/duncanwa/akari/data_test")
HDF5_PATH = Path("/mn/stornext/d16/cmbco/bp/gustavbe/master/dirbe_hdf5_files/")
CIO_PATH = Path('/mn/stornext/d16/cmbco/ola/akari/TODs')
CIO_PATH = Path('/mn/stornext/d16/cmbco/ola/akari/TODs/manually_extracted')

# system constants
N_PROC = multiprocessing.cpu_count()

ROTATOR = hp.Rotator(coord=["E", "G"])
YDAYS = np.concatenate([np.arange(89345, 89366), np.arange(90001, 90265)])
START_TIME = Time("1981-01-01", format="isot", scale="utc")

# CIO constants
N_CIO_FILES = 285
BAD_DATA_SENTINEL = -16375
TSCAL = 2e-15
SAMP_RATE = 1 / 26
#25.28
#16.86

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
        tods = []
        pixels = []
        flags = []
        time_starts = []
        time_stops = []
        sat_pos_start = []
        sat_pos_stop = []
        earth_pos_start = []
        earth_pos_stop = []
        for yday in yday_data:
            if f'AKARI_{band}' in yday.tods.keys():
                tods.append(yday.tods[f"AKARI_{band}"])
                pixels.append(yday.pixels[f"AKARI_{band}"])
                flags.append(yday.flags[f"AKARI_{band}"])
                time_starts.append(yday.time_start)
                time_stops.append(yday.time_stop)
                sat_pos_start.append(yday.earth_pos_start)
                sat_pos_stop.append(yday.earth_pos_stop)
                earth_pos_start.append(yday.earth_pos_start)
                earth_pos_stop.append(yday.earth_pos_stop)
        output_dict[f'AKARI_{band}'] =  CIO(tods=tods,
             pixels=pixels,
             flags=flags,
             time_start=time_starts,
             time_stop=time_stops,
             sat_pos_start=sat_pos_start,
             sat_pos_stop=sat_pos_stop,
             earth_pos_start=earth_pos_start,
             earth_pos_stop=earth_pos_stop)
    return output_dict


def get_yday_data(
        tod: NDArray, lon: NDArray, lat: NDArray, flags: NDArray, nside_out: int, planet_time_delta: timedelta, color_corr: bool, 
) -> list[YdayData]:


    if len(tod) == 100:
        # Processing SW pixels

        bands = ['WIDE-S']*60 + ['N60']*40
        inds = np.concatenate((np.arange(60), np.arange(40))) + 1
        data_ind = np.arange(100)
    else:
        # Processing LW pixels

        bands = ['WIDE-L']*45 + ['N160']*30
        inds = np.concatenate((np.arange(45), np.arange(30))) + 1
        data_ind = np.arange(75)

    with multiprocessing.Pool(processes=N_PROC) as pool:
        proc_chunks = [
            pool.apply_async(
                get_yday_cio_data,
                args=(tod[di], lon[di], lat[di], flags[di], band, ind, nside_out, color_corr),
            )
            for band, ind, di  in zip(bands, inds, data_ind)
        ]
        return [result.get() for result in proc_chunks if result]


def get_yday_cio_data(
    tod: NDArray,
    lon: NDArray,
    lat: NDArray,
    flag: NDArray,
    band: str,
    ind: int,
    nside_out: int,
    #planet_interps: dict[str, dict[str, interp1d]],
    color_corr: bool,
) -> YdayData:
    """Function which extracts and reorders the CIO data from one day CIO file."""


    #yday = YDAYS[file_number]

    t = 53826

    # Convert time to MJD
    t = (START_TIME + TimeDelta(t, format="sec", scale="tai")).mjd
    sat_pos_start, earth_pos_start = 0,0
    sat_pos_stop, earth_pos_stop = 0,0


    # Extract tods, and modify pointing vectors per detector according to beam data
    tods: dict[str, np.ndarray] = {}
    pixels: dict[str, np.ndarray] = {}
    flags: dict[str, np.ndarray] = {}
    band_label = f'AKARI_{band}_{ind:02}'



    c = SkyCoord(ra=lon*u.deg, dec=lat*u.deg, frame='icrs')
    pixels[band_label] = hp.ang2pix(nside_out, c.galactic.l.value, c.galactic.b.value, lonlat=True)

    tods[band_label] = tod




    flags[band_label] = flag

    return YdayData(
        tods, 
        pixels, 
        flags, 
        time_start=t, 
        time_stop=t, 
        sat_pos_start=sat_pos_start, 
        sat_pos_stop=sat_pos_stop, 
        earth_pos_start=earth_pos_start,
        earth_pos_stop=earth_pos_stop,
    )


def padd_array_gaps(splits: list[np.ndarray], padding: list[np.ndarray]) -> np.ndarray:
    return np.concatenate([np.append(s, p) for s, p in zip(splits, padding)])





def write_band(
        comm_tod: TODLoader, cio: CIO, filename: str, band: str, ndet: int, nside_out: int, n_pids: int
        , pid_0: int, raw: bool) -> None:
    COMMON_GROUP = "/common"
    HUFFMAN_COMPRESSION = ["huffman", {"dictNum": 1}]
    HUFFMAN_COMPRESSION_TOD = ["huffman", {"dictNum": 2}]

    det_str = ""
    for i in range(ndet):
        det_str = det_str + f'AKARI_{band}_{i+1:02}'
        if i < ndet - 1:
            det_str = det_str + ', '
    comm_tod.init_file(freq=filename, od="", mode="w")



    for v in cio.keys():
        det1 = v
        break

    if ('N60' in band) or ('WIDE-S' in band):
        comm_tod.add_field(COMMON_GROUP + "/fsamp", 25.28)
    else:
        comm_tod.add_field(COMMON_GROUP + "/fsamp", 16.86)
    comm_tod.add_field(COMMON_GROUP + "/nside", [nside_out])
    comm_tod.add_field(COMMON_GROUP + "/det", np.string_(det_str + ","))

    comm_tod.add_field(COMMON_GROUP + "/polang", [0]*ndet)
    comm_tod.add_attribute(COMMON_GROUP + "/polang", "index", det_str)

    comm_tod.add_field(COMMON_GROUP + "/mbang", [0]*ndet)
    comm_tod.add_attribute(COMMON_GROUP + "/mbang", "index", det_str)


    for pid in range(n_pids):
        pid_label = f"{pid+pid_0:06}"
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
            det_lab = f'AKARI_{band}_{i+1:02}'
            pid_det_group = f"{pid_label}/{det_lab}"
            comm_tod.add_field(pid_det_group + "/flag", cio[det_lab].flags[pid], HUFFMAN_COMPRESSION)

            if raw:
                comm_tod.add_field(pid_det_group + "/ztod", cio[det_lab].tods[pid], HUFFMAN_COMPRESSION_TOD)
            else:
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
    chip: str,
    pid_0: int,
    raw: bool,
    overwrite: bool = False,
) -> None:
    manager = multiprocessing.Manager()
    pool = multiprocessing.Pool(processes=N_PROC)

    multiprocessor_manager_dicts = {}
    filenames = {}
    for band in akari_utils.BANDS:
        # name = f"DIRBE_{band:02}_nside{nside_out:03}_zodi_only"
        name = f"AKARI_{band}_n{nside_out}_v{version:02}"
        multiprocessor_manager_dicts[name] = manager.dict()
        filenames[band] = name

    comm_tod = TODLoader(
        out_path, "", version, multiprocessor_manager_dicts, overwrite=overwrite
    )

    n_pids = 0
    for key, cio in zip(cios.keys(), cios.values()):
        n_pids = len(cio.time_start)
        break
    if n_pids == 0:
        raise ValueError("No CIOs found")

    # Currently, writing commander h5 files uses a pretty unoptimal interface which is difficult to
    # parallelize without corrupting the written files. This is the go to way to concurrently write files
    # which is pretty ugly.

    if chip == 'LW':
        NDETS = akari_utils.NDETS[2:]
        BANDS = akari_utils.BANDS[2:]
    elif chip == 'SW':
        NDETS = akari_utils.NDETS[:2]
        BANDS = akari_utils.BANDS[:2]
    x = [[]]
    for ndet, band in zip(akari_utils.NDETS, akari_utils.BANDS):
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
                n_pids,
                pid_0,
                raw)))

    for res1 in x:
        for res in res1:
            res.get()

    pool.close()
    pool.join()

    comm_tod.make_filelists()


# SW fsamp 24 Hz, LW 16 Hz
# Normal sampling, they reset either 0.5 sec, 1 sec, or 2 sec, to avoid detector
# saturation.
# Differential sampling/CDS, they only take a sample every two points in high
# intensity regions.

# We should make sure that the data are split based on the reset timing.

# Once a minute, the calibration lamp is flashed for less than one second.

def main() -> None:
    time_delta = timedelta(hours=1)

    pid_0 = 1
    files = range(1)
    #nside_out = 512   # 2**9
    #nside_out = 8192  # 2**13
    nside_out = 2048  # 2**11

    start_time = time.perf_counter()
    color_corr = False
    version = 3

    raw = True


    from glob import glob
    #fnames = glob(CIO_PATH/f'flux/SW/FIS_SW_??????????????_flux.pkl')

    fnames = list(CIO_PATH.glob('flux/SW/FIS_SW_*_flux.pkl'))
    fnames.sort()

    print(fnames)


    #fnames = fnames[:2]
    #fnames = fnames[13:15]

    times = []
    for f in fnames:
        times.append(str(f).split('FIS_SW_')[1][:14])


    print(f"{'Writing AKARI h5 files':=^50}")
    print(f"{version=}, {nside_out=}")
    print(f"reading and processing cios for {len(files)} ydays...")
    pid_now = 0
    yday_data = []
    for t in times:
        pid_now += 1
        for chip in ['LW', 'SW']:
            if raw:
                datfile = CIO_PATH/f'ADU/{chip}/FIS_{chip}_{t}_adu.pkl'
            else:
                datfile = CIO_PATH/f'flux/{chip}/FIS_{chip}_{t}_flux.pkl'
            latfile = CIO_PATH/f'lat/{chip}/FIS_{chip}_{t}_gb_lat.pkl'
            lonfile = CIO_PATH/f'lon/{chip}/FIS_{chip}_{t}_gb_lon.pkl'
            with open(datfile, 'rb') as f:
                tods = pickle.load(f)
                if raw:
                    tods = np.array(tods).astype(int)
            with open(latfile, 'rb') as f:
                lats = pickle.load(f)
            with open(lonfile, 'rb') as f:
                lons = pickle.load(f)


            flag_ind = 0

            bad_frames = CIO_PATH/f'frame_flag/bad_frame/{chip}/FIS_{chip}_{t}_flame.pkl'
            with open(bad_frames, 'rb') as f:
                flag = np.array(pickle.load(f)).astype(int)
                flag_ind += 1

            flag_tot = np.zeros(flag.shape, dtype=int)

            pixel_frames = []
            flag_list = ['dead', 'no_diff', 'reset', 'rstanom', 'saturate', 'gpgl_tail']

            flag_list += [f'gpgl_type{i}' for i in range(1,5)]
            flag_list += [f'mtgl_type{i}' for i in range(1,5)]


            for flag in flag_list:
                pixel_frames.append(CIO_PATH/f'pixel_flag/{flag}/{chip}/FIS_{chip}_{t}_gb_{flag}.pkl')

            for pf in pixel_frames:
                with open(pf, 'rb') as f:
                    flag = np.array(pickle.load(f)).astype(int)
                    try:
                        flag_tot[flag != 0] += 2**(flag_ind)
                    except IndexError:
                        print(flag.shape, pf, 'not included because shape is wrong')
                    flag_ind += 1


            status_list = ['calalon', 'calason', 'calbon', 'shtop', 'sinalon', 'sinason']
            status_frames = []
            for sl in status_list:
                status_frames.append(CIO_PATH/f'status/{sl}/{chip}/FIS_{chip}_{t}_gb_{sl}.pkl')

            for status, stat in zip(status_list, status_frames):
                with open(stat, 'rb') as f:
                    flag = np.array(pickle.load(f)).astype(int)
                    if 'shtop' == status:
                        flag = 1-flag
                    try:
                        flag_tot[flag != 0] += 2**(flag_ind)
                    except IndexError:
                        print(flag.shape, pf, 'not included because shape is wrong')
                    flag_ind += 1

            stat = CIO_PATH/f'no_peri_corr/{chip}/FIS_{chip}_{t}_gb_no_peri_corr.pkl'
            with open(stat, 'rb') as f:
                flag = np.array(pickle.load(f)).astype(int)
                try:
                    flag_tot[flag != 0] += 2**(flag_ind)
                except IndexError:
                    print(flag.shape, pf, 'not included because shape is wrong')
                flag_ind += 1

            





            yday_data += get_yday_data(
                tods, lons, lats, flag_tot, nside_out=nside_out, planet_time_delta=time_delta, color_corr=color_corr,
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
            pid_0=pid_0,
            chip=chip,
            raw=raw,
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
