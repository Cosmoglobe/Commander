import pickle
import h5py
from tqdm import tqdm
from functools import cache

import numpy as np
import healpy
import glob
import astropy.units as u
from pathlib import Path
from astropy.io import fits
import matplotlib.pyplot as plt 
import shutil
import os

from iras_native_tod_reader import BANDS, BAND_DETS
from commander_npz2hdf import run_multiple_band_write_external



def akari_beam():
    akari = "/mn/stornext/d23/cmbco/globe/akari"
    beam = f"{akari}/beam"
    ff = glob.glob(f"{beam}/*.fits")
    ff = sorted(ff)

def akari_bp():
    akari = "/mn/stornext/d23/cmbco/globe/akari"
    bp = h5py.File(f"{akari}/bp/AKARI_instrument_v05.h5")
    akari_bands = ["065", "090", "140", "160"]
    fwhm_list = [37, 39, 58, 61]
    fwhm = {k:v*u.arcsec for k,v in zip(akari_bands, fwhm_list)}

    labels = {
        "065": "N60",
        "090": "WIDE-S",
        "140": "WIDE-L",
        "160": "N160",
    }
    labels
    colors = ["b", "purple", "r", "brown"]
    ms = [10,5,10,5]
    for i, band in enumerate(akari_bands):
        # print(f"fwhm / {band} = {fwhm[band].to(u.arcmin)}")
        bandy = bp[f"AKARI_{band}"]["bandpass"][:]
        bandx = bp[f"AKARI_{band}"]["bandpassx"][:]
        bandx = (bandx * u.GHz).to(u.micron, equivalencies=u.spectral())
        plt.plot(bandx, bandy, label=labels[band], c=colors[i])
        plt.plot(bandx[0], bandy[0], 'o', ms=ms[i], c=colors[i])
        plt.plot(bandx[-1], bandy[-1], 's', ms=ms[i], c=colors[i])


        print(f"{band}: {len(bandx)=}. bandx = [{bandx[:5]}, {bandx[-5:]}]")
        print(f"{band}: {len(bandy)=}. bandy = [{bandy[:5]}, {bandy[-5:]}]")

        # print(f"{band}: {len(bandx)} {len(bandx[::-1])}")
    plt.legend()
    plt.xlabel("micron")
    plt.xlim(30, 240)
    plt.ylim(-0.2, 1.4)
    plt.grid()
    plt.show()

    # print(f"{bsandx}")
    # print()s
    # print(f"{bsandy}")

    # plt.plot()




def bad_chunks():
    # KEPT FOR REFERENCE ONLY.
    # Currently not being used.
    bdets_60 = {
        0:  8 , # chunk 5: len(data) = 8933 
        1:  9 , # chunk 5: len(data) = 8933 
        2:  10, # chunk 5: len(data) = 8933 
        3:  11, # chunk 5: len(data) = 8933 
        4:  12, #### chunk 5: len(data) = 8932 / 12
        5:  13, #### chunk 5: len(data) = 8932 / 13
        6:  14, #### chunk 5: len(data) = 8932 / 14
        7:  15, #### chunk 5: len(data) = 8932 / 15
        8:  31, # chunk 5: len(data) = 8933
        9:  32, # chunk 5: len(data) = 8933
        10: 33, # chunk 5: len(data) = 8933
        11: 34, # chunk 5: len(data) = 8933
        12: 35, #### chunk 5: len(data) = 8932 / 35
        13: 37, #### chunk 5: len(data) = 8932 / 37
        14: 38  #### chunk 5: len(data) = 8932 / 38
    }

    bdets_100 = {
        0:  1 , # chunk 5: len(data) = 4467
        1:  2 , # chunk 5: len(data) = 4467
        2:  3 , # chunk 5: len(data) = 4467
        3:  4 , # chunk 5: len(data) = 4467
        4:  5 , # chunk 5: len(data) = 4467
        5:  6 , # chunk 5: len(data) = 4467
        6:  7 , # chunk 5: len(data) = 4466
        7:  55, # chunk 5: len(data) = 4466
        8:  56, # chunk 5: len(data) = 4467
        9:  57, # chunk 5: len(data) = 4467
        10: 58, # chunk 5: len(data) = 4467
        11: 59, # chunk 5: len(data) = 4467
        12: 60, # chunk 5: len(data) = 4467
        13: 61, # chunk 5: len(data) = 4466
        14: 62  # chunk 5: len(data) = 4466
    } 






def test_tod_files(version, nproc=256, nside=512):
    """
    Make tod files for a given version with external script.
    Both h5 files and filelists will be stored in:
        **/globe/iras/tod/vetle_deplated_IPAC/v{version}/n{nside}
    These are NOT read by commander by default.
    Commander reads tod files from filelists located at
        filelist_dir = **/globe/iras/tod/{band}/fileslist_{band}.txt
    
    copy_new_filelists: 
        True:   Copy new h5 filelists to filelist_dir
        False:  Make new h5 files, but don't replace filelists. 
    """

    test_outpath = f"/mn/stornext/d5/data/vetleav/IRAS/tod_data/hdf_commander/test{version}"

    bands = BANDS
    bands = BANDS[:1]

    band_dets_copy = BAND_DETS.copy()
    use_every = 1
    # band_dets_list = BAND_DETS[bands[0]]
    # band_dets_list.remove("IRAS_012-25")

    # band_dets = {band: band_dets_copy[band][::use_every] for band in bands}
    if version == 6:
        band_dets = {band: band_dets_copy[band][3:] for band in bands}
    if version == 7:
        band_dets = {band: band_dets_copy[band][2:-1] for band in bands}
    if version == 8:
        band_dets = {band: band_dets_copy[band][1:-2] for band in bands}
    if version == 9:
        band_dets = {band: band_dets_copy[band][0:-3] for band in bands}
    if version == 10:
        band_dets = {band: band_dets_copy[band] for band in bands}
        band_dets[bands[0]].pop(1)
    if version == 11:
        band_dets = {band: band_dets_copy[band][1:3] for band in bands}

    # band_dets = {band: band_dets_copy[band] for band in bands}
    # band_dets = ban
    print(f"Storing TODs:\n    ver: {version}\n    dst: {test_outpath}")

    print("BANDS: ",*band_dets.items(), sep="\n")
    # print(f"{len(band_dets[bands[0]])}")
    # print()
    # print(os.cpu_count())
    # print(2**9)

    # return
    # exit()

    run_multiple_band_write_external(
        outpath         = test_outpath,
        version         = version,
        nside           = nside,
        overwrite       = True,
        num_processes   = nproc,#os.cpu_count(),
        bands           = bands,
        band_dets       = band_dets,
        num_files_per_segment = 100
    )
    
test_tod_files(version=11, nproc=256)
exit()
test_tod_files(version=6, nproc=256)
test_tod_files(version=7, nproc=256)
test_tod_files(version=8, nproc=os.cpu_count())
test_tod_files(version=9, nproc=os.cpu_count())
