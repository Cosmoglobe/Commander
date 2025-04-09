"""
Script to compare the current FIRAS maps to the Planck official maps in the corresponding frequencies.
"""
import os

import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits

path = "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/fits_files/"
planck_path = "planck/"
firas_path = "maps/frequency_maps/rh_ss/galactic/"

save_path = "/mn/stornext/d16/www_cmb/aimartin/firas/maps/comparison/"
if not os.path.exists(save_path):
    os.makedirs(save_path)

bands = ["070", "100", "143", "217", "353", "545", "857"]
firas_bands = ["0068", "0095", "0149", "0217", "0353", "0544", "0857"]

for i, band in enumerate(bands):
    # if int(band) > 639:
    #     firas_path = "maps/frequency_maps/rh_ss/galactic/"
    planck_filename = "planck_" + band + "_smoothed.fits"
    firas_filename = f"{firas_bands[i]}_offset0.5_nside32.fits"

    planck_map = fits.open(path + planck_path + planck_filename)[0].data
    firas_map = fits.open(path + firas_path + firas_filename)[0].data

    # check if firas_map is nonetype
    if firas_map is None:
        print(firas_filename + " is None")
        continue
    diff = planck_map - firas_map
    hp.mollview(diff, title="Planck - FIRAS " + band + " GHz", nest=False, xsize=2000, min=-20, max=20)
    hp.graticule()
    plt.savefig(save_path + band + ".png")