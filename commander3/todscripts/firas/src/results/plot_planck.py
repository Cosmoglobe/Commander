"""
Script to check that the smoothed Planck maps are correctly generated.
"""

import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
grandparent = os.path.dirname(parent)
sys.path.append(grandparent)
import globals as g

planck = False
firas = True

user = os.environ["USER"]

if planck:
    path = f"/mn/stornext/u3/{user}/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/fits_files/planck/"
    bands = ["070", "100", "143", "217", "353", "545", "857"]

    save_path = f"{g.SAVE_PATH}maps/planck/"

    for band in bands:
        map = fits.open(path + "planck_" + band + "_smoothed.fits")[0].data
        hp.mollview(map, title="Planck " + band + " smoothed", nest=False, xsize=2000, min=0, max=200)
        hp.graticule()
        plt.savefig(save_path + "planck_" + band + "_smoothed.png")

if firas:
    path = f"/mn/stornext/u3/{user}/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/fits_files/"
    bands = ["0068", "0095", "0149", "0217", "0353", "0544", "0857"]

    save_path = f"{g.SAVE_PATH}maps/planck_equivalent/"

    channel = "rh_ss"

    for band in bands:
        # if int(band) > 639:
        #     channel = "rh_ss"
        map = fits.open(path + f"maps/frequency_maps/{channel}/galactic/" + band + "_offset0.5_nside32.fits")[0].data
        hp.mollview(map, title=band + " GHz", nest=False, xsize=2000, min=0, max=200)
        hp.graticule()
        plt.savefig(save_path + band + ".png")
