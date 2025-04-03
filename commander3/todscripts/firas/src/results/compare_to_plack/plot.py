"""
Script to check that the smoothed Planck maps are correctly generated.
"""

import os

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

planck = False
firas = True

user = os.environ["USER"]

if planck:
    path = f"/mn/stornext/u3/{user}/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/fits_files/planck/"
    bands = ["070", "100", "143", "217", "353", "545", "857"]

    save_path = f"/mn/stornext/d16/www_cmb/{user}/firas/maps/planck/"
    # check if path exists and create it if it doesn't
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    for band in bands:
        map = fits.open(path + "planck_" + band + "_smoothed.fits")[0].data
        hp.mollview(map, title="Planck " + band + " smoothed", nest=False, xsize=2000, min=0, max=200)
        hp.graticule()
        plt.savefig(save_path + "planck_" + band + "_smoothed.png")

if firas:
    path = f"/mn/stornext/u3/{user}/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/fits_files/"
    bands = ["0068", "0095", "0149", "0217", "0353", "0544", "0857"]

    save_path = f"/mn/stornext/d16/www_cmb/{user}/firas/maps/planck_equivalent/"
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    channel = "rh_ss"

    for band in bands:
        # if int(band) > 639:
        #     channel = "rh_ss"
        map = fits.open(path + f"maps/frequency_maps/{channel}/galactic/" + band + "_offset0.5_nside32.fits")[0].data
        hp.mollview(map, title=band + " GHz", nest=False, xsize=2000, min=0, max=200)
        hp.graticule()
        plt.savefig(save_path + band + ".png")
