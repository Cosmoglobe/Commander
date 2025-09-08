"""
Script to compare the current FIRAS maps to the Planck official maps in the corresponding frequencies.
"""

import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
grandparent = os.path.dirname(parent)
sys.path.append(grandparent)
import globals as g

path = "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/fits_files/"
planck_path = "planck/"
firas_path = "maps/frequency_maps/"

save_path = f"{g.SAVE_PATH}maps/comparison/"

bands = [70, 100, 143, 217, 353, 545, 857]
firas_bands = [68, 95, 149, 217, 353, 544, 857]

for i, band in enumerate(bands):
    # if int(band) > 639:
    #     firas_path = "maps/frequency_maps/rh_ss/galactic/"
    planck_filename = f"planck_{band:03d}_smoothed.fits"
    if band < 600:
        channel_mode = "ll_ss/"
    else:
        channel_mode = "rh_ss/"
    firas_filename = f"{channel_mode}galactic/{firas_bands[i]:04d}_nside32.fits"

    planck_map = fits.open(path + planck_path + planck_filename)[0].data
    firas_map = fits.open(path + firas_path + firas_filename)[0].data

    # check if firas_map is nonetype
    if firas_map is None:
        print(firas_filename + " is None")
        continue
    diff = planck_map - firas_map
    hp.mollview(
        diff,
        title=f"Planck - FIRAS {band:04d} GHz",
        nest=False,
        xsize=2000,
        min=-20,
        max=20,
        cmap="RdBu_r",
    )
    hp.graticule()
    plt.savefig(f"{save_path}diff_{band}.png")

    ratio = planck_map / firas_map
    hp.mollview(
        ratio,
        title=f"Planck / FIRAS {band:04d} GHz",
        nest=False,
        xsize=2000,
        min=0,
        max=2,
        cmap="RdBu_r",
    )
    hp.graticule()
    plt.savefig(f"{save_path}ratio_{band}.png")
