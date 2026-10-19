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
import os
from multiprocessing import Pool
from glob import glob
import multiprocessing as mp

N_PROC = multiprocessing.cpu_count()
ROTATOR = hp.Rotator(coord=["E", "G"])
YDAYS = np.concatenate([np.arange(89345, 89366), np.arange(90001, 90265)])
START_TIME = Time("1981-01-01", format="isot", scale="utc")

N_CIO_FILES = 285
BAD_DATA_SENTINEL = -16375
TSCAL = 2e-15
SAMP_RATE = 1 / 26

SAMP_RATE_DAYS = SAMP_RATE / (24 * 3600)

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
    for band in akari_utils.DETECTORS: #detまで入る
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
                # print(f"yday.tods[AKARI_{band}]=", yday.tods[f"AKARI_{band}"]) #-> 値入ってる
                tods.append(yday.tods[f"AKARI_{band}"]) #1detの[1時間分のデータ, 1時間分のデータ,...]
                pixels.append(yday.pixels[f"AKARI_{band}"])
                flags.append(yday.flags[f"AKARI_{band}"])
                # print("yday.time_start", yday.time_start) #-> 値入ってる(全部DIRBEのtime_startだけど)
                time_starts.append(yday.time_start)
                time_stops.append(yday.time_stop)
                sat_pos_start.append(yday.earth_pos_start)
                sat_pos_stop.append(yday.earth_pos_stop)
                earth_pos_start.append(yday.earth_pos_start)
                earth_pos_stop.append(yday.earth_pos_stop)
        # if 'N60' in f'AKARI_{band}':
        #     print(f"AKARI_{band} tods=", tods)

        # if ('N60' in f'AKARI_{band}') or ('WIDE-S' in f'AKARI_{band}'):
        #     if tods == []:
        #         print(f"AKARI_{band} tods=", tods)
        # else:
        #     if tods != []:
        #         print(f"AKARI_{band} tods=", tods) 
        # -> ok
        
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

    if len(tod) == 100:  # Processing SW pixels
        bands = ['N60']*40 + ['WIDE-S']*60
        inds = np.concatenate((np.arange(40), np.arange(60))) + 1
        data_ind = np.arange(100)
    else:  # Processing LW pixels
        bands = ['WIDE-L']*45 + ['N160']*30
        inds = np.concatenate((np.arange(45), np.arange(30))) + 1
        data_ind = np.arange(75)

    results = []
    for band, ind, di in zip(bands, inds, data_ind):
        result = get_yday_cio_data(tod[di], lon[di], lat[di], flags[di], band, ind, nside_out, color_corr)
        results.append(result)

    return results

def get_yday_cio_data(
    tod: NDArray,
    lon: NDArray,
    lat: NDArray,
    flag: NDArray,
    band: str,
    ind: int,
    nside_out: int,
    color_corr: bool,
) -> YdayData:

    t = 53826

    # Convert time to MJD
    t = (START_TIME + TimeDelta(t, format="sec", scale="tai")).mjd
    # print(t)
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
        if band in v:
            # print(v)
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

    # print("n_pids "+ str(n_pids))
    for pid in range(n_pids):
        # print("pid "+str(pid))
        pid_label = f"{pid+pid_0:06}"
        # print("pid_label "+str(pid_label))
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
                print(pid_det_group + "/ztod", cio[det_lab].tods[pid])
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
    BAND: str,
    pid_0: int,
    raw: bool,
    outputname: str,
    overwrite: bool = False,
) -> None:
    multiprocessor_manager_dicts = {}
    filenames = {}
    for band in akari_utils.BANDS:
        name = outputname
        multiprocessor_manager_dicts[name] = {}  # manager.dict()の代わりに通常の辞書を使用
        filenames[band] = name

    comm_tod = TODLoader(
        out_path, "", version, multiprocessor_manager_dicts, overwrite=overwrite
    )

    n_pids = 1

    # n_pids = 0
    # # print(cios.keys())
    # # print(len(cios.values()))
    # for key, cio in zip(cios.keys(), cios.values()):
    #     n_pids = len(cio.time_start)
    #     # print("cio.time_start "+str(cio.time_start))
    #     break
    # if n_pids == 0:
    #     raise ValueError("No CIOs found")

    # if chip == 'LW':
    #     NDETS = akari_utils.NDETS[2:]
    #     BANDS = akari_utils.BANDS[2:]
    # elif chip == 'SW':
    #     NDETS = akari_utils.NDETS[:2]
    #     BANDS = akari_utils.BANDS[:2]

    # for ndet, band in zip(akari_utils.NDETS, akari_utils.BANDS):
    #     # print()
    #     write_band(
    #         comm_tod,
    #         cios,
    #         filenames[band],
    #         band,
    #         ndet,
    #         nside_out,
    #         n_pids,
    #         pid_0,
    #         raw
    #     )
    band = BAND
    if band == 'N60':
        # print(band, akari_utils.BANDS[0])
        ndet = akari_utils.NDETS[0]
    elif band == 'WIDE-S':
        print(band, akari_utils.BANDS[1])
        ndet = akari_utils.NDETS[1]
    elif band == 'WIDE-L':
        # print(band, akari_utils.BANDS[2])
        ndet = akari_utils.NDETS[2]
    elif band == 'N160':
        # print(band, akari_utils.BANDS[3])
        ndet = akari_utils.NDETS[3]
    write_band(
        comm_tod,
        cios,
        filenames[band],
        band,
        ndet,
        nside_out,
        n_pids,
        pid_0,
        raw
    )
    
    filelist = []
    for freq in multiprocessor_manager_dicts.keys():
        for buf in multiprocessor_manager_dicts[freq].values():
            filelist.append(buf)

    comm_tod.make_filelists()

    return filelist[0]


def process_t(t, raw, nside_out, time_delta, color_corr, version, pid_0, AKARI_DATA_PATH, BAND, CIO_PATH):
    h = t[:-4]
    # print("h: "+str(h))
    filelists = []
    yday_data = []
    
    # for chip in ['LW', 'SW']:
    if (BAND == 'N60') or (BAND == 'WIDE-S'):
        chips = ['SW']
    elif (BAND == 'WIDE-L') or (BAND == 'N160'):
        chips = ['LW']
    for chip in chips:
        # print(chip)

        # if raw:  # det, flux, lat, lon
        #     datfile = CIO_PATH/f'det/all/{h}/FIS_{chip}_{t}_det.pkl'
        # else:
        datfile = CIO_PATH/f'flux/all/{h}/FIS_{chip}_{t}_flux_re.pkl'
        latfile = CIO_PATH/f'lat/all/{h}/FIS_{chip}_{t}_gb_lat_gads.pkl'
        lonfile = CIO_PATH/f'lon/all/{h}/FIS_{chip}_{t}_gb_lon_gads.pkl'

        with open(datfile, 'rb') as f:
            tods = pickle.load(f)
            if raw:
                tods = np.array(tods).astype(float).astype(int)
        with open(latfile, 'rb') as f:
            lats = pickle.load(f)
        with open(lonfile, 'rb') as f:
            lons = pickle.load(f)

        flag_ind = 0  # bad_frame
        bad_frames = CIO_PATH/f'flag/all/{h}/FIS_{chip}_{t}_bad_frame.pkl'
        with open(bad_frames, 'rb') as f:
            flag = np.array(pickle.load(f)).astype(int)
            # flag_ind += 1
            flag_ind += 0

        flag_tot = np.zeros(flag.shape, dtype=int)

        pixel_frames = []  # pixelのflag
        flag_list = ['bad','dead','saturate','reset_re','rstanom','mtgl_tail','flutter','no_rp_corr','blank']
        # flag_list += [f'gpgl_type{i}' for i in range(1,5)]
        flag_list += [f'mtgl_type{i}' for i in range(1,5)]

        for flag in flag_list:
            pixel_frames.append(CIO_PATH/f'pix_flag/all/{h}/FIS_{chip}_{t}_{flag}.pkl')

        for pf in pixel_frames:
            with open(pf, 'rb') as f:
                flag = np.array(pickle.load(f)).astype(int)
                try:
                    flag_tot[flag != 0] += 2**(flag_ind)
                except IndexError:
                    print(flag.shape, pf, 'not included because shape is wrong')
                flag_ind += 1
                
        status_list = ['calalon_re', 'calason_re', 'shtop']  # status
        status_frames = [CIO_PATH/f'status/all/{h}/FIS_{chip}_{t}_{sl}.pkl' for sl in status_list]

        for status, stat in zip(status_list, status_frames):
            with open(stat, 'rb') as f:
                flag = np.array(pickle.load(f)).astype(int)
                if 'shtop' == status:
                    flag = 1-flag
                try:
                    flag_tot[flag != 0] += 2**(flag_ind)
                except IndexError:
                    print(flag.shape, 'not included because shape is wrong')
                flag_ind += 1

        yday_data += get_yday_data(
            tods, lons, lats, flag_tot, nside_out=nside_out, planet_time_delta=time_delta, color_corr=color_corr,
        )

    cios = get_cios(yday_data)

    # print(cios)

    cio_time = time.perf_counter()
    # print(f"Finished reading: {h}")
    # print("done")
    # print("writing cios to h5 files...")

    file_exp = write_to_commander_tods(
        cios,
        nside_out=nside_out,
        version=version,
        out_path=AKARI_DATA_PATH,
        pid_0=pid_0,
        BAND=BAND,
        raw=raw,
        overwrite=True,
        outputname=str(t)+'_flux_reflag',
    )

    filelists.append(file_exp)

    print(f"Finished writing cios to h5: ", chip, BAND, h)

    return filelists

# 並列処理を実行する関数
def parallel_processing(times, raw, nside_out, time_delta, color_corr, version, pid_0, AKARI_DATA_PATH, BAND, CIO_PATH):
    num_processes = mp.cpu_count()  # 利用するCPUコア数を設定
    with mp.Pool(processes=80) as pool:
        results = pool.starmap(process_t, [(t, raw, nside_out, time_delta, color_corr, version, pid_0, AKARI_DATA_PATH, BAND, CIO_PATH) for t in times])
        # 作成したhdf
    
    # 結果をまとめる
    filelists = []
    for result in results:
        filelists.extend(result)
    
    # print(filelists)
    return filelists

def multi(BAND, AKARI_DATA_PATH, CIO_PATH):
    time_delta = timedelta(hours=1)

    pid_0 = 1
    files = range(1)
    nside_out = 2048  # 2**11

    start_time = time.perf_counter()
    color_corr = False
    version = 3

    raw = True

    dir_path = Path("lat/all")
    files_dir = [f for f in os.listdir(CIO_PATH / dir_path) if os.path.isdir(CIO_PATH / dir_path / f)]

    fnames = []
    for path in files_dir:
        fnames.append(list((CIO_PATH / dir_path / path).glob('FIS_SW_*.pkl')))

    fnames.sort()

    times = []
    for f in fnames:
        times.append(str(f).split('FIS_SW_')[1][:14])
    # times = [times[2616]]
    print(times)

    start_time = time.perf_counter()
    filelists = parallel_processing(times, raw, nside_out, time_delta, color_corr, version, pid_0, AKARI_DATA_PATH, BAND, CIO_PATH)

    # # path_w = "/mn/stornext/d23/cmbco/cg/AKARI/ryosukem2/data/akari/filelist.txt"
    # path_w = f"/mn/stornext/d23/cmbco/cg/AKARI/tamakim3/data/data_{BAND}/filelist.txt"
    # num_lines = len(filelists)  # 書き込む行数を計算

    # # print(filelists)

    # with open(path_w, mode='w') as f:
    #     f.write(f"{num_lines}\n")  # 最初の行に行数を書き込む
    #     f.writelines(filelists)    # filelistsの内容を書き込む

    # total_time = time.perf_counter() - start_time
    # print(f"Total processing time: {(total_time / 60):2.2f} minutes")

if __name__ == '__main__':
    BAND = 'N160'
    #AKARI_DATA_PATH = Path("/mn/stornext/d23/cmbco/cg/AKARI/ryosukem/data/akari")
    AKARI_DATA_PATH = Path(f"/mn/stornext/d23/cmbco/cg/AKARI/tamakim3/data/data_{BAND}")
    CIO_PATH = Path('/mn/stornext/d23/cmbco/cg/AKARI/akari_TSD_pkl')

    multi(BAND, AKARI_DATA_PATH, CIO_PATH)
